.. This is the top-level walkthrough of the AstrOptical package.

Background
==========

This package arose from the need to process several different data types, both
observational and simulated, within the visible astronomy field. Originally,
the tools found here were a conglomeration of scripts I used to import and
process data. However, they quickly grew to be too large and scattered.
Therefore, I packaged them up and published them hoping they'd save time for
someone else in the future.

Overview
========

AstrOptical comes packaged with several different data processing and data I/O
tools. Many of the tools here are built off of
`Astropy <http://www.astropy.org/>`_ tools. For example, the cosmicray module
simply wraps the `Astroscrappy <https://github.com/astropy/astroscrappy>`_
package to keep the user from having to perform redundant tasks of file I/O.

On the other hand, the tool also adds methods to read in data from tools such
as `SExtractor <https://www.astromatic.net/software/sextractor>`_ and
`Starburst99 <http://www.stsci.edu/science/starburst99/docs/default.htm>`_.
Again, I could not find a good toolset out there to read these data in Python
so I decided to write my own.
