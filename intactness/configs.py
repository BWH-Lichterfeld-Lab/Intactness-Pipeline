"""Set up basic running parameters

Note: See configparser https://docs.python.org/3/library/configparser.html
"""

import os
import logging
import datetime
from configparser import ConfigParser, ExtendedInterpolation

# pylint: disable=C0103
# Invalid constant name
logger = logging.getLogger('pipe')


def configs(filename):
    """Parse configurations
    """
    logger.info('Checking parameters')

    cfg = ConfigParser(interpolation=ExtendedInterpolation())
    cfg.read(filename)

    if cfg['Main']['email'] == '':
        print('Please specify your email address')
        exit(1)

    folder = cfg['Main']['path_out']
    if os.path.exists(folder):
        msg = '{} already existed'.format(folder)
        logger.info(msg)
    else:
        os.makedirs(folder)

    # Create a new log file each time the program is started
    date_tag = datetime.datetime.now().strftime("%Y-%b-%d_%H-%M-%S")
    file_log = "{}/run_{}.log".format(cfg['Main']['path_out'], date_tag)

    # Logging format: access time, logger name and message
    fmtr = logging.Formatter('%(asctime)s\t\t%(name)s\t\t%(message)s')

    fh = logging.FileHandler(file_log)
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmtr)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(fmtr)
    logger.addHandler(ch)

    return cfg
