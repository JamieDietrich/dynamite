# Python Process Runner
# Author: Jerry Dietrich
# 2022 June 8

from multiprocessing.pool import Pool
import multiprocessing
import threading
import importlib
import traceback
import inspect
import socket
import sys
import os
    
EXECUTE_lock = multiprocessing.Semaphore()
PPR_instance = None

def init(lock):
    global pipe_lock
    pipe_lock = lock

class NoDaemonProcess(multiprocessing.Process):
    @property
    def daemon(self):
        return False
    @daemon.setter
    def daemon(self, value):
        pass

class NoDaemonContext(type(multiprocessing.get_context())):
    Process = NoDaemonProcess

class NestablePool(Pool):
    def __init__(self, *args, **kwargs):
        kwargs['context'] = NoDaemonContext()
        super(NestablePool, self).__init__(*args, **kwargs)

class PPR(object):
    def __init__(self, instance_name, processes = None, remote_terminate = True):
        global PPR_instance
        if instance_name != None:
            if PPR_instance == None:
                PPR_instance = self
                self.pool_dict = {}
                self.parent_conn, self.child_conn = multiprocessing.Pipe()
                thread = threading.Thread(target=self.create_pipes, args=())
                thread.daemon = True
                thread.start()
                if remote_terminate:
                    thread = threading.Thread(target=self.terminate_PPR, args=(True, ))
                    thread.daemon = True
                    thread.start()
                self.cpu_count = multiprocessing.cpu_count() if processes == None else processes
                self.pool_dict = {'main': NestablePool(self.cpu_count, initializer = init, initargs = (multiprocessing.Lock(),))}
            else:
                self.remove_non_main_pools()
            if isinstance(instance_name, str):
                PPR_instance.instance = (getattr(importlib.import_module(instance_name), instance_name), False)
            elif isinstance(instance_name, tuple):
                if isinstance(instance_name[0], str) and isinstance(instance_name[1], str):
                    PPR_instance.instance = (getattr(importlib.import_module(instance_name[0]), instance_name[1]), False)
                else:
                    PPR_instance.instance = (getattr(instance_name[0], instance_name[1]), False) if instance_name[1] != None else (instance_name[0], True)
            else:
                print("Invalid Instance Name")    
            
    def __new__(self, *args, **kwargs):
        return super().__new__(self) if PPR_instance == None else PPR_instance
    
    def terminate_PPR(self, remote_terminate = False):
        global PPR_instance
        if remote_terminate:
            try:
                kill_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                kill_socket.bind(("127.0.0.1", 65000))
                kill_socket.listen()
                kill_socket.accept()[0].close()
                kill_socket.close()
                if self.child_conn == None: return
            except Exception as e:
                print(e)
                print("Remote Terminate not possible because socket is already in use or has an error") 
                return
        try:
            self.parent_conn.close()
            self.child_conn.close()
            self.child_conn = None
            PPR_instance = None
            for pool in self.pool_dict.values():
                pool.close() 
                pool.terminate()
            try:
                with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                    s.setdefaulttimeout(5)
                    s.connect(("127.0.0.1", 65000))
            except: pass
        except: pass
        if remote_terminate: os._exit(0)
    
    def remove_non_main_pools(self):
        pool = self.pool_dict['main']
        for p in self.pool_dict.keys():
            if p != 'main':
                self.pool_dict[p].close() 
                self.pool_dict[p].terminate()
        self.pool_dict = {'main': pool}
                       
    def create_pipes(self):
        try:
            while self.child_conn != None:
                pool = self.child_conn.recv()
                if pool not in self.pool_dict.keys():
                    self.pool_dict[pool] = multiprocessing.Pool(self.cpu_count)
                parent_conn, child_conn = multiprocessing.Pipe()
                thread = threading.Thread(target=self.create_processes_jobs, args=(child_conn, self.pool_dict[pool]))
                thread.daemon = True
                thread.start()
                self.child_conn.send(parent_conn)
        except Exception as e:
            print(e)
        
    def create_processes_jobs(self, c_conn, pool):
        with c_conn:
            c_conn.send(self.create_processes(*c_conn.recv(), pool = pool))
    
    def run_functions(self, params):
        return self.run_functions_ppr(params, False)
     
    def run_functions_ppr(self, params, ppr = True):
        ((instance, func), args) = params
        try:
            if ppr:
                return instance(*args,self) if not inspect.isclass(instance) else (getattr(instance(), func)(*args,self))
            else:
                return instance(*args) if not inspect.isclass(instance) else (getattr(instance(), func)(*args))
        except:
            print(traceback.format_exc())
            try:socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect(("127.0.0.1", 65000))
            except: pass
        finally:
            sys.stdout.flush()
            sys.stderr.flush()

    def create_processes(self, function = None, args = None, size = None, callback = None, pool = None, ppr = False):
        def get_func(func):
            return (getattr(self.instance[0], func), self.instance[0]) if self.instance[1] else (self.instance[0], func)

        if function == None: return None
        if not hasattr(self, 'pool_dict') or self.pool_dict == {}:
            with pipe_lock:
                self.parent_conn.send(function[0] if isinstance(function, list) else function)
                p_conn = self.parent_conn.recv()
            with p_conn:
                p_conn.send((function, args, size, callback))
                return p_conn.recv()
        if pool == None: pool = self.pool_dict['main']
        if size != None and size != 0 and not isinstance(function, list):
            results = pool.map(self.run_functions_ppr if ppr else self.run_functions,[(get_func(function), t) for t in self.process_begin_size_tuple(args, (0, size))])
        elif isinstance(args, list):
            results = pool.map(self.run_functions_ppr if ppr else self.run_functions,[(get_func(function[i]) if isinstance(function, list) else get_func(function), args[i]) for i in range(len(args))])
        elif isinstance(function, list):
            results = pool.map(self.run_functions_ppr if ppr else self.run_functions,[(get_func(function[i]), args) for i in range(len(function))])
        else:
            if ppr:
                results = getattr(self.instance[0]() if not self.instance[1] else self.instance[0], function)(*args,self)
            else:
                results = getattr(self.instance[0]() if not self.instance[1] else self.instance[0], function)(*args)
        if results == None or len(results) == 0: return [] if callback == None else callback([])
        return results if callback == None else callback(results)
               
    def execute_in_parallel(self, function, *args):
        try:
            EXECUTE_lock.acquire()
            #return getattr(self.instance[0], function)(self.instance[0], *(args if isinstance(args, tuple) else tuple(args)))
            return getattr(self.instance[0], function)(*(args if isinstance(args, tuple) else tuple(args)))
        except:
            print(traceback.format_exc())
        finally:
            EXECUTE_lock.release()
                
    def process_begin_size_tuple(self, arg, params):
        (begin, size) = params
        args = []
        start = begin
        end = size
        if size < 0:
            for s in range(size*-1):
                args.append((arg + (s,)))
        elif (size > self.cpu_count):
            for cpu in range(self.cpu_count):
                if (cpu == 0):
                    start = begin
                    end = (int)(size/self.cpu_count) + begin
                elif (cpu == self.cpu_count -1): 
                    start = ((int)(size/self.cpu_count) * cpu) + begin
                    end = size + begin
                else:
                    start = ((int)(size/self.cpu_count) * cpu) + begin
                    end = ((int)(size/self.cpu_count) * (cpu + 1)) + begin
                args.append((arg + (start,) + (end,)))
        else:
            args.append((arg + (start,) + (size,)))
        return args
            
    def __getstate__(self):
        self_dict = self.__dict__.copy()
        try: del self_dict['pool_dict']
        except: pass
        return self_dict

    def __setstate__(self, state):
        self.__dict__.update(state)
        
if __name__ == '__main__':
    PPR(None)
