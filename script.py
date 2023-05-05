import os
import sys
from irgsctool import GenerateIRGSC as GC
from irgsctool import ValidateIRGSC as VG


gc = VG(66.75, 15.86)
gc.validate(validate=True)