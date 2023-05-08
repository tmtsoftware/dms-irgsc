import os
import sys
from irgsctool import GenerateIRGSC as GC
from irgsctool import ValidateIRGSC as VG


gc = VG(334.27, 0.38)
gc.validate(validate=True)