"""
Exceptions for the primertool package.
"""


class PrimertoolError(Exception):
    """Base class for exceptions in the primertool module."""
    pass


class PrimertoolGenomeError(PrimertoolError):
    """Exception raised for errors in the genome."""
    pass


class PrimertoolInputError(PrimertoolError):
    """Exception raised for errors in the input."""
    pass


class PrimertoolMutalyzerError(PrimertoolError):
    """Exception raised for errors in the Mutalyzer API."""
    pass


class PrimertoolExonLengthError(PrimertoolError):
    """Exception raised when the exon length exceeds the max insert size."""
    pass


class PrimertoolNoPrimerFoundError(PrimertoolError):
    """Exception raised when no primer is found."""
    pass


class PrimertoolIntronicPositionError(PrimertoolError):
    """Exception raised when the given position is intronic."""
    pass


class InSilicoPCRError(object):
    """Base class for exceptions in insilicopcr module."""
    pass


class InSilicoPCRRequestError(InSilicoPCRError):
    """Exception raised request errors to the InSilicoPCR API."""
    pass
