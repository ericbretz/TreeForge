from typing import Optional, Dict, Any
from pathlib import Path
from core.utils.printout import PrintOut


class BaseStage:
    def __init__(
        self,
        log: int,
        hc: str,
        bc: str,
        threads: Optional[int] = None,
        subprocess_dir: Optional[Path] = None,
        shared_printClass: Optional[PrintOut] = None
    ):
        self.log            = log
        self.hc             = hc
        self.bc             = bc
        self.threads        = threads
        self.subprocess_dir = subprocess_dir

        if shared_printClass is not None:
            self.printClass = shared_printClass
        else:
            self.printClass = PrintOut(log, hc, bc)
        self.printout       = self.printClass.printout
        
        self.return_dict: Dict[str, Any] = {}
    
    def run(self):
        raise NotImplementedError("subclasses must implement this")
