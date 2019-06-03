import logging
import os
import sys


def init_logger(log_path, log_file, print_log=True, level=logging.INFO):
    if not os.path.isdir(log_path):
        os.makedirs(log_path)

    fileHandler = logging.FileHandler("{0}/{1}.log".format(log_path, log_file))

    handlers = [fileHandler]

    if print_log:
        consoleHandler = logging.StreamHandler(sys.stdout)
        handlers.append(consoleHandler)

    logging.basicConfig(
        level=level,
        format=
        "%(asctime)s [%(process)d] [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s",
        handlers=handlers)
