"""
This module contains the logger configuration for the primertool package.
"""
import logging

logger = logging.getLogger(__package__)


class CustomFormatter(logging.Formatter):
    """Logging colored formatter, adapted from https://stackoverflow.com/a/56944256/3638629

    Attributes:
        fmt (str): The format string to use for the log messages
        FORMATS (dict): A dictionary mapping log levels to their respective format strings
    """

    def __init__(self, fmt: str):
        """
        Initialize the formatter with a format string.

        Args:
            fmt (str): The format string to use for the log messages
        """
        super().__init__()

        grey = '\x1b[38;21m'
        blue = '\x1b[38;5;39m'
        yellow = '\x1b[38;5;226m'
        red = '\x1b[38;5;196m'
        bold_red = '\x1b[31;1m'
        reset = '\x1b[0m'

        self.fmt = fmt
        self.FORMATS = {
            logging.DEBUG: grey + self.fmt + reset,
            logging.INFO: blue + self.fmt + reset,
            logging.WARNING: yellow + self.fmt + reset,
            logging.ERROR: red + self.fmt + reset,
            logging.CRITICAL: bold_red + self.fmt + reset
        }

    def format(self, record: logging.LogRecord) -> str:
        """Format the log record according to its level.

        Args:
            record (logging.LogRecord): The log record to format
        Returns:
            str: The formatted log message
        """
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


def init_logger(fmt: str = '[%(levelname)s] %(message)s', level: int = logging.DEBUG,
                save_to: str = None) -> logging.Logger:
    """Initialize the logger for the primertool package.

    Args:
        fmt (str): The format string to use for the log messages
        level (int): The logging level to use
        save_to (str): The file to save the logs to (optional)
    Returns:
        logging.Logger: The configured logger
    """
    logger = logging.getLogger(__package__)
    if save_to:
        logging.basicConfig(filename=save_to)
    logger.setLevel(level)

    # Create stdout handler for logging to the console (logs all five levels)
    stdout_handler = logging.StreamHandler()
    stdout_handler.setLevel(level)
    stdout_handler.setFormatter(CustomFormatter(fmt))

    # Add the stdout handler to the logger
    logger.addHandler(stdout_handler)

    return logger
