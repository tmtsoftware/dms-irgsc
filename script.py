import os
import sys
from irgsctool import GenerateIRGSC as GC
from irgsctool import ValidateIRGSC as VG


gc = VG(35, -3.5)
gc.validate(validate=True)