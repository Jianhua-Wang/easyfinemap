"""Main module."""

import logging


class EasyFinemap(object):
    """Main class."""

    def __init__(self):
        """Initialize."""
        self.logger = logging.getLogger(__name__)

    def run(self):
        """Run the program."""
        self.logger.info("Running the program")
