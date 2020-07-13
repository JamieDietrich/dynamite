from itertools import chain
import multiprocessing
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
        if self.instance[1]: return (getattr(self.instance[0], func), None)
        return (self.instance[0], func)
    
    def create_processes(self, function, args, size = None, return_sorted = True, block = True):
        processes, results = [], []
        if (size != None):
            results = getattr(self.pool,"map" if block else "apply_async")(run_functions,[(self.get_func(function), t) for t in self.process_begin_size_tuple(args, (0, size))],)
            if not block: results = []
        else:
            if isinstance(args, list):
                for i in range(len(args)):
                    processes.append(self.pool.apply_async(run_functions, [(self.get_func(function[i] if isinstance(function, list) else function), args[i]),]))
                if block:
                    for i in range(len(processes)):
                        results.append(processes[i].get())
            else:
                results = getattr(self.pool,"map" if block else "apply_async")(run_functions, [(self.get_func(function), args)])
                if not block: results = []
        #if not return_sorted: return ", ".join(repr(e) for e in results)
        return results
        results.sort(key = lambda d: list(d.keys())[0])
        results = [list(chain([it for its in i for it in its])) for i in zip(*chain.from_iterable(map(lambda d: list(d.values()),results)))]
        return ", ".join(repr(e) for e in results[0]) if len(results) == 1 else ", ".join(repr(e) for e in results)
        
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
    
