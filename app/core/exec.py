import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from queue import Queue
from collections import Counter
from threading import Lock
from config.__init__ import echo


def pool(tasks, max_workers, iteration=''):
    # tasks = ['cmd-1', 'cmd-2']
    # results = [{'res': 'stdout', 'code': 0}, {'res': 'stdout', 'err': 'stderr', 'code': 1}]
    task_queue = Queue()
    for task in tasks:
        task_queue.put(task)

    stats = Counter(total=len(tasks), completed=0)
    stats_lock = Lock()
    results = []

    def info():
        # remain = stats['total'] - stats['completed']
        spc = (f" ({iteration})" if iteration else "") + (" " * 20)
        echo(f"➜ Runs: {stats['completed']}/{stats['total']}{spc}\r", 36)

    info()

    def task_wrapper(task):
        try:
            # log += f"Task: {task}\n"
            res = subprocess.run(
                task, shell=True, capture_output=True, text=True)
            with stats_lock:
                stats['completed'] += 1
            return task, res.returncode, res.stdout, res.stderr
        except Exception as e:
            with stats_lock:
                stats['completed'] += 1
            raise e

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        while not task_queue.empty() and len(futures) < max_workers:
            task = task_queue.get()
            futures.append(executor.submit(task_wrapper, task))

        while futures:
            for future in as_completed(futures):
                futures.remove(future)
                cmd, code, stdout, stderr = future.result()
                results.append((cmd, code, stdout, stderr,))
                # print("Done", cmd, code)
                if code != 0:
                    # log += f"Error: {stderr}\n"
                    print(f"-> {cmd}: {code}")
                    print(f"<- Error: {stderr}")
                if not task_queue.empty():
                    next_task = task_queue.get()
                    futures.append(executor.submit(task_wrapper, next_task))

            with stats_lock:
                info()

    echo("\r  " + (" " * 100) + "\r")
    return results
