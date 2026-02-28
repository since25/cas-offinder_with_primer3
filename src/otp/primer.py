import primer3
from typing import Dict, Any, Optional

class PrimerDesigner:
    def __init__(self, global_args: Optional[Dict[str, Any]] = None):
        """
        Initialize the PrimerDesigner with default or user-provided primer3 arguments.
        """
        self.global_args = global_args or {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 5,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 8,
            'PRIMER_MAX_SELF_END': 3,
            'PRIMER_PAIR_MAX_COMPL_ANY': 8,
            'PRIMER_PAIR_MAX_COMPL_END': 3,
            'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]]
        }
    
    def design(self, flank_seq: str, target_start_in_flank0: int, target_len: int, num_return: int = 1, size_range: list = [150, 250]) -> dict:
        """
        Design primers for a given sequence to encompass the target region.
        
        Variables use 0-based indexing. Target region is [target_start_in_flank0, target_start_in_flank0 + target_len).
        """
        # Ensure a 20bp buffer around the target so primers don't overlap or sit exactly on the edge
        buffer = 20
        start = max(0, target_start_in_flank0 - buffer)
        length = target_len + (2 * buffer)
        
        seq_args = {
            'SEQUENCE_ID': 'offtarget_flank',
            'SEQUENCE_TEMPLATE': flank_seq,
            'SEQUENCE_TARGET': [start, length]
        }
        
        # Override user runtime arguments
        args = self.global_args.copy()
        args['PRIMER_NUM_RETURN'] = num_return
        args['PRIMER_PRODUCT_SIZE_RANGE'] = [size_range]
        
        def run_primer3(current_args):
            try:
                return primer3.bindings.design_primers(seq_args, current_args)
            except Exception as e:
                return {"error": str(e), "success": False, "PRIMER_PAIR_NUM_RETURNED": 0}
        
        results = run_primer3(args)
        
        # If no primers returned, try relaxing the parameters
        if results.get('PRIMER_PAIR_NUM_RETURNED', 0) == 0:
            relaxed_args = args.copy()
            relaxed_args.update({
                'PRIMER_MIN_TM': 53.0,
                'PRIMER_MAX_TM': 67.0,
                'PRIMER_MIN_GC': 15.0,
                'PRIMER_MAX_GC': 85.0,
                'PRIMER_MAX_POLY_X': 6,
                'PRIMER_MAX_NS_ACCEPTED': 1,
                'PRIMER_PAIR_MAX_COMPL_ANY': 10,
                'PRIMER_PAIR_MAX_COMPL_END': 5
            })
            results = run_primer3(relaxed_args)
            
        if results.get('PRIMER_PAIR_NUM_RETURNED', 0) == 0:
            return {"error": results.get("error", "No primers returned even with relaxed parameters"), "success": False}
            
        # Parse the best (first) primer pair
        left_pos_len = results.get('PRIMER_LEFT_0', [-1, -1])
        right_pos_len = results.get('PRIMER_RIGHT_0', [-1, -1])
        
        left_seq = results.get('PRIMER_LEFT_0_SEQUENCE', '')
        right_seq = results.get('PRIMER_RIGHT_0_SEQUENCE', '')
        
        left_tm = results.get('PRIMER_LEFT_0_TM', 0.0)
        right_tm = results.get('PRIMER_RIGHT_0_TM', 0.0)
        
        penalty = results.get('PRIMER_PAIR_0_PENALTY', 0.0)
        amplicon_size = results.get('PRIMER_PAIR_0_PRODUCT_SIZE', 0)
        
        covers = False
        if left_pos_len[0] != -1 and right_pos_len[0] != -1:
            # right primer position returned by primer3 is the 3' end. So its coordinate is reversed.
            # Its 5' position is right_pos_len[0] - right_pos_len[1] + 1
            # For coverage checking:
            # left primer 3' end = left_pos_len[0] + left_pos_len[1] - 1
            # right primer 3' end = right_pos_len[0]
            # left covers target if left primer end <= target_start
            left_end = left_pos_len[0] + left_pos_len[1]
            right_start = right_pos_len[0] - right_pos_len[1] + 1
            if left_end <= target_start_in_flank0 and right_start >= target_start_in_flank0 + target_len:
                covers = True
        
        return {
            "success": True,
            "primer_left_pos_in_flank0": left_pos_len[0],
            "primer_left_len": left_pos_len[1],
            "primer_right_pos_in_flank0": right_pos_len[0],
            "primer_right_len": right_pos_len[1],
            "primer_left_seq": left_seq,
            "primer_right_seq": right_seq,
            "primer_left_tm": left_tm,
            "primer_right_tm": right_tm,
            "amplicon_size": amplicon_size,
            "covers_offtarget": covers,
            "primer_pair_penalty": penalty
        }
