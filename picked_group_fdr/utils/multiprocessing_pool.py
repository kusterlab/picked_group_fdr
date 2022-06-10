from __future__ import print_function

import sys
import time
import signal
import warnings
import logging
from multiprocessing import Pool
from logging.handlers import QueueHandler, QueueListener
import multiprocessing.pool
import traceback

logger = logging.getLogger(__name__)


# https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic
class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)


# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
# https://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic#54304172
class NestablePool(multiprocessing.pool.Pool):
    def Process(self, *args, **kwds):
        proc = super(NestablePool, self).Process(*args, **kwds)
        proc.__class__ = NoDaemonProcess

        return proc


class JobPool:
    def __init__(self, processes=1, warningFilter="default", queue=None):
        self.warningFilter = warningFilter
        
        # In the GUI, SIMSI-Transfer runs in a child process to allow user
        # interaction with the GUI itself. In this case, we need a 
        # NoDaemonProcess to allow this child process to spawn children,
        # which is implemented in the NestablePool class.
        # We also need to pass the logger to these grandchildren processes
        # using the multiprocessing.Queue and QueueListener classes.
        if not queue and multiprocessing.current_process().name != 'MainProcess':
            queue = multiprocessing.Queue()
            queue_listener = QueueListener(queue, logger)
            queue_listener.start()
        self.pool = NestablePool(processes, worker_init, initargs=(self.warningFilter, queue))
        
        self.results = []

    def applyAsync(self, f, fargs, *args, **kwargs):
        r = self.pool.apply_async(f, fargs, *args, **kwargs)
        self.results.append(r)

    def checkPool(self, printProgressEvery=-1):
        try:
            outputs = list()
            for res in self.results:
                outputs.append(res.get(timeout=10000))  # 10000 seconds = ~3 hours
                if printProgressEvery > 0 and len(outputs) % printProgressEvery == 0:
                    logger.info(
                        f' {len(outputs)} / {len(self.results)} {"%.2f" % (float(len(outputs)) / len(self.results) * 100)}%')
            self.pool.close()
            self.pool.join()
            return outputs
        except (KeyboardInterrupt, SystemExit):
            logger.error("Caught KeyboardInterrupt, terminating workers")
            self.pool.terminate()
            self.pool.join()
            sys.exit()
        except Exception as e:
            logger.error("Caught Unknown exception, terminating workers")
            logger.error(traceback.print_exc())
            logger.error(e)
            self.pool.terminate()
            self.pool.join()
            sys.exit()

    def stopPool(self):
        self.pool.terminate()
        self.pool.join()


def worker_init(warningFilter, queue=None):
    if queue:
        queueHandler = QueueHandler(queue)
        logger = logging.getLogger()
        logger.setLevel(logging.INFO)
        logger.addHandler(queueHandler)

    # set warningFilter for the child processes
    warnings.simplefilter(warningFilter)

    # causes child processes to ignore SIGINT signal and lets main process handle
    # interrupts instead (https://noswap.com/blog/python-multiprocessing-keyboardinterrupt)
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def addOne(i):
    return i + 1


def unitTest():
    pool = JobPool(4)
    for i in range(20):
        pool.applyAsync(addOne, [i])
    results = pool.checkPool()
    print(results)
