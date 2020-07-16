import multiprocessing
import collections
import importlib
import traceback
import os

def run_functions(params):
    ((instance, func), args) = params
    try:
        return (getattr(instance(), func)(*args) if func != None else instance(*args))
    except KeyboardInterrupt:
        os._exit(1);
    except Exception:
        print(traceback.format_exc())
            
class Client:
    def __init__(self, instance_name, processes = None):
        if (instance_name != None):
            self.cpu_count = multiprocessing.cpu_count() if processes == None else processes
            self.pool = multiprocessing.Pool(processes)
            if isinstance(instance_name, str):
                self.instance = (getattr(importlib.import_module(instance_name), instance_name), False)
            elif isinstance(instance_name, tuple):
                self.instance = (getattr(instance_name[0], instance_name[1]), False) if instance_name[1] != None else (instance_name[0], True)
            else:
                print("Invalid Instance Name")
                
    def get_func(self, func):
        return (getattr(self.instance[0], func), None) if self.instance[1] else (self.instance[0], func)
    
    def create_processes(self, function, args, size = None, callback = None, async = False):
        if (size != None):
            rslt = getattr(self.pool,"apply_async" if async else "map")(run_functions,[(self.get_func(function), t) for t in self.process_begin_size_tuple(args, (0, size))],)
        else:
            if isinstance(args, list):
                processes, rslt = [], []
                for i in range(len(args)):
                    processes.append(self.pool.apply_async(run_functions, [(self.get_func(function[i] if isinstance(function, list) else function), args[i]),]))
                if not async:
                    for r in range(len(processes)): rslt.append(processes[r].get())
            else:
                rslt = getattr(self.pool,"apply_async" if async else "map")(run_functions, [(self.get_func(function), args)])
        if async: return {} if callback == None else callback({})
        results = collections.OrderedDict((key,d[key]) for d in sorted(rslt, key = lambda d: list(d.keys())[0]) for key in d)
        return results if callback == None else callback(results)
    
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
        del self_dict['pool']
        return self_dict

    def __setstate__(self, state):
        self.__dict__.update(state)
        
if __name__ == '__main__':
    Client(None)
