"""
logger.py

Logging module for cavar

Can be used by calling:
    from logger import logMessage as logger

Example:
    logger("warning", f"File {file }is missing")
"""
import sys
import logging

# Configure the logger
_logger = logging.getLogger(__name__)
_logger.setLevel(logging.DEBUG)

# Create console handler
_console_handler = logging.StreamHandler(stream = sys.stderr)
_console_handler.setLevel(logging.DEBUG)

# Set format
_formatter = logging.Formatter(fmt = "[%(asctime)s] %(levelname)s %(message)s", datefmt = "%Y-%m-%d %H:%M:%S")

# Add formatting to console handler
_console_handler.setFormatter(_formatter)

# Add console handler to logger
if not _logger.handlers:
    _logger.addHandler(_console_handler)

def log_message(level: str, message: str) -> None:
    """
    Portable function to expose logging to the
    main script. Required level and message can
    be passed to it for debug. Always writes to
    STDERR.

    Levels: debug, info, warning, error, critical
    """
    level = level.lower()
    log_func = getattr(_logger, level, None)

    if callable(log_func):
        log_func(message)
    else:
        _logger.warning("Unknown level '%s', message: %s", level, message)
