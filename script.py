import os
import sys
from irgsctool import GenerateIRGSC as GC
from irgsctool import ValidateIRGSC as VG


gc = VG(269.93, -13.48)
gc.validate(validate=True)