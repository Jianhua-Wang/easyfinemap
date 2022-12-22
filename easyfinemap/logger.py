"""Define custom logger."""

import logging

from rich.logging import RichHandler

logging.basicConfig(
    level=logging.NOTSET,
    format="%(name)s - %(levelname)s - %(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
)

logger = logging.getLogger(__name__)
