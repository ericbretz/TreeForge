import os
import yaml
from pathlib             import Path
from core.utils.printout import PrintOut

class ConfigManager:    
    def __init__(self, highlight_color, background_color):
        self.highlight_color  = highlight_color
        self.background_color = background_color
        self.config_dir       = Path.home() / '.treeforge'
        self.printClass       = PrintOut('', self.highlight_color, self.background_color)
        self.printout         = self.printClass.printout
        
    def get_defaults_dict(self):
        return {
            'input_dir'                : os.getcwd(),
            'iter'                     : 5,
            'threads'                  : 2,
            'clutter'                  : False,
            'output_dir'               : '',
            'blast_evalue'             : 10.0,
            'blast_max_targets'        : 1000,
            'mcl_hit_frac_cutoff'      : 0.3,
            'mcl_minimum_taxa'         : 10,
            'mcl_inflation'            : 1.4,
            'mcl_perfect_identity'     : 100.0,
            'mcl_coverage_threshold'   : 0.9,
            'mcl_min_seq_length'       : 300,
            'mafft_maxiter'            : 1000,
            'mafft_pxclsq_threshold'   : 0.1,
            'mafft_thread_divisor'     : 4,
            'tree_relative_cutoff'     : 0.2,
            'tree_absolute_cutoff'     : 0.3,
            'tree_branch_cutoff'       : 0.02,
            'tree_mask_paralogs'       : 'n',
            'tree_outlier_ratio'       : 20.0,
            'tree_max_trim_iterations' : 10,
            'tree_min_subtree_taxa'    : 4,
            'tree_min_leaves'          : 4,
            'prune_orthocutoff'        : 20,
            'prune_relative_cutoff'    : 0.2,
            'prune_absolute_cutoff'    : 0.3,
            'prune_outlier_ratio'      : 20.0,
            'prune_max_trim_iterations': 10,
            'prune_min_tree_leaves'    : 3,
            'prank_seqtype'            : 'aa',
            'prank_pxclsq_threshold'   : 0.3,
            'prank_bootstrap'          : 1000,
            'super_bootstrap'          : 1000,
            'output_super_matrix'      : False,
        }
    
    def get_parameter_types(self):
        return {
            'input_dir'                : 'string',
            'iter'                     : 'integer',
            'threads'                  : 'integer',
            'clutter'                  : 'boolean',
            'output_dir'               : 'string',
            'blast_evalue'             : 'float',
            'blast_max_targets'        : 'integer',
            'mcl_hit_frac_cutoff'      : 'float',
            'mcl_minimum_taxa'         : 'integer',
            'mcl_inflation'            : 'float',
            'mcl_perfect_identity'     : 'float',
            'mcl_coverage_threshold'   : 'float',
            'mcl_min_seq_length'       : 'integer',
            'mafft_maxiter'            : 'integer',
            'mafft_pxclsq_threshold'   : 'float',
            'mafft_thread_divisor'     : 'integer',
            'tree_relative_cutoff'     : 'float',
            'tree_absolute_cutoff'     : 'float',
            'tree_branch_cutoff'       : 'float',
            'tree_mask_paralogs'       : 'string',
            'tree_outlier_ratio'       : 'float',
            'tree_max_trim_iterations' : 'integer',
            'tree_min_subtree_taxa'    : 'integer',
            'tree_min_leaves'          : 'integer',
            'prune_orthocutoff'        : 'integer',
            'prune_relative_cutoff'    : 'float',
            'prune_absolute_cutoff'    : 'float',
            'prune_outlier_ratio'      : 'float',
            'prune_max_trim_iterations': 'integer',
            'prune_min_tree_leaves'    : 'integer',
            'prank_seqtype'            : 'string',
            'prank_pxclsq_threshold'   : 'float',
            'prank_bootstrap'          : 'integer',
            'super_bootstrap'          : 'integer',
            'output_super_matrix'      : 'boolean',
        }
    
    def load_config(self, config_path):
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def _extract_values(self, config):
        values = {}
        
        for section_name, section in config.items():
            if section_name == 'treeforge_config':
                continue
                
            if isinstance(section, dict):
                if 'description' in section:
                    for param_name, param_data in section.items():
                        if param_name != 'description' and isinstance(param_data, dict):
                            if 'value' in param_data:
                                values[param_name] = param_data['value']
                else:
                    for key, value in section.items():
                        if key != 'treeforge_config':
                            values[key] = value
        return values
    
    def save_config(self, config, output_path):
        self.config_dir.mkdir(exist_ok=True)
        
        with open(output_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False, indent=2)
    
    def create_config(self, output_path):
        config = {
            **self.get_defaults_dict()
        }
        self.save_config(config, output_path)
        out_name = str(output_path)
        out_name = out_name if len(out_name) < 32 else '...' + out_name[-29:]
        self.printout('info', f"Config created at: {out_name}")
    
    def validate_config(self, config):
        parameter_types = self.get_parameter_types()
        defaults        = self.get_defaults_dict()
        
        for key, value in config.items():
            if key not in defaults:
                self.printout('error', f"Unknown config parameter: {key}")
                continue
            
            expected_type = parameter_types.get(key, 'any')
            if expected_type == 'integer' and not isinstance(value, int):
                self.printout('error', f"{key} must be an integer, got {type(value).__name__}")
                return False
            elif expected_type == 'float' and not isinstance(value, (int, float)):
                self.printout('error', f"{key} must be a number, got {type(value).__name__}")
                return False
            elif expected_type == 'boolean' and not isinstance(value, bool):
                self.printout('error', f"{key} must be a boolean, got {type(value).__name__}")
                return False
            elif expected_type == 'string' and not isinstance(value, str):
                self.printout('error', f"{key} must be a string, got {type(value).__name__}")
                return False
        
        return True 