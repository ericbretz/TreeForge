import os
import sys
import ete3
import shutil
import subprocess
from pathlib import Path
from core.stages.base_stage import BaseStage
from core.utils.sublogger import run_stage_subprocess, run_logged_subprocess
from core.utils.bes import BES
from core.utils.concat import SequenceConcatenator

class Astral(BaseStage):
    def __init__(self,
                 dir_prank,
                 dir_super,
                 dir_treeforge,
                 renamed_map,
                 output_super_matrix,
                 super_bootstrap,
                 threads,
                 log,
                 hcolor,
                 bcolor,
                 coalescent_tree_path,
                 molecular_tree_path,
                 supermatrix_tree_path,
                 gene_trees_dir,
                 id_to_stem=None,
                 subprocess_dir=None,
                 shared_printClass=None):
        super().__init__(log, hcolor, bcolor, threads, subprocess_dir, shared_printClass)
        self.dir_prank           = Path(dir_prank)
        self.dir_super           = Path(dir_super)
        self.dir_treeforge       = Path(dir_treeforge)
        self.renamed_map         = renamed_map
        self.id_to_stem          = id_to_stem if id_to_stem is not None else {}
        self.output_super_matrix = output_super_matrix
        self.super_bootstrap     = super_bootstrap

        self.cln_files           = ' '.join([str(self.dir_prank / f) for f in os.listdir(self.dir_prank) if f.endswith('.treefile')])
        self.alignment_files     = ' '.join([str(self.dir_prank / f) for f in os.listdir(self.dir_prank) if f.endswith('-cln')])
        self.concat_tree         = str(self.dir_super / 'concat.tre')
        self.finaltree           = str(coalescent_tree_path)
        self.supermatrix         = str(self.dir_super / 'SuperMatrix.treefile')
        self.genetree            = str(supermatrix_tree_path)
        self.dir_gene_trees      = Path(gene_trees_dir)
        self.molecular_tree_path = molecular_tree_path
        self.s_matrix_dir        = str(self.dir_super / 'super.matrix')
        self.s_model_dir         = str(self.dir_super / 'super.model')
        self.return_dict         = {'concat_tree'    : self.concat_tree,
                                    'gene_trees'     : self.dir_gene_trees,
                                    'finaltree'      : self.finaltree,
                                    'num_trees'      : len(self.cln_files.split()),
                                    'supermatrix_taxa': None,
                                    'supermatrix_bp'  : None,
                                    'bes_tree_count'  : None}

    def run(self):
        """
        Run the final stage of the pipeline.
        """
        self.file_check()
        self.tree_concat()
        self.astral()
        self.bes()
        if self.output_super_matrix:
            self.concatenate_alignments()
            self.iqtree()

        return self.return_dict
    
    def file_check(self):
        """
        Check if the gene trees files exist.
        """
        if len(self.cln_files.split()) == 0:
            self.printout('error', 'No gene trees found')
            self.printout('error', 'Try adjusting parameters')
            sys.exit(1)
        if self.output_super_matrix and len(self.alignment_files.split()) == 0:
            self.printout('error', 'No alignment files found for supermatrix')
            self.printout('error', 'Try adjusting parameters')
            sys.exit(1)
    
    def tree_concat(self):
        """
        Concatenate gene trees and rename leaves to original file names.
        Saves both a directory of gene trees and a concatenated tree.
        """
        self.printout('metric', 'Concatenating gene trees')
        self.dir_gene_trees.mkdir(parents=True, exist_ok=True)
        processed_files = 0
        unmapped_count  = 0
        
        with open(self.concat_tree, 'w') as outfile:
            for file in self.cln_files.split():
                try:
                    tree = ete3.Tree(file)
                    tree_gene = tree.copy(method='deepcopy')
                    
                    for node in tree.traverse():
                        if node.is_leaf():
                            lookup_key = node.name.replace('_', '@', 1)  # Replace only first underscore
                            if lookup_key in self.renamed_map:
                                _, source_file = self.renamed_map[lookup_key]
                                node.name = source_file.split('.')[0]
                            else:
                                node.name = node.name.split('.')[0]
                                unmapped_count += 1
                    
                    for node in tree_gene.traverse():
                        if node.is_leaf():
                            lookup_key = node.name.replace('_', '@', 1)  # Replace only first underscore
                            if lookup_key in self.renamed_map:
                                original_name, _ = self.renamed_map[lookup_key]
                                node.name = original_name
                            else:
                                node.name = node.name
                    
                    tree_gene.write(format=1, outfile=str(self.dir_gene_trees / (file.split('/')[-1].split('_')[0] + '.tre')))
                    outfile.write(tree.write(format=1) + '\n')
                    processed_files += 1
                except Exception as e:
                    self.printout('error', f'Failed to process tree file {file}: {str(e)}')
                    continue
        
        # self.printout('metric', f'Successfully processed {processed_files} tree files')
        
        if unmapped_count > 0:
            self.printout('warning', f'{unmapped_count} gene IDs could not be mapped to species names')
            self.printout('warning', 'This may cause ASTRAL to produce an incorrect species tree')

    def astral(self):
        """
        Run ASTRAL.
        ASTRAL creates a coalescent species tree from the concatenated gene trees.
        """
        # self.printout('metric', 'Astral')
        if Path(self.concat_tree).exists():
            with open(self.concat_tree, 'r') as f:
                content = f.read().strip()
                if len(content) == 0:
                    self.printout('error', 'Concatenated tree file is empty')
                    sys.exit(1)
        else:
            self.printout('error', f'Concatenated tree file not found: {self.concat_tree}')
            sys.exit(1)

        self.printout('metric', 'Running ASTRAL on concatenated gene trees')
        cmd = f'astral -i {self.concat_tree} -o {self.finaltree}'
        
        run_stage_subprocess(cmd, 'ASTRAL', self.subprocess_dir, self.printout)

    def concatenate_alignments(self):
        """
        Concatenate sequences to create a supermatrix.
        """
        self.printout('metric', 'Creating Supermatrix')
        file_list_path = self.dir_super / 'alignment_file_list.txt'
        with open(file_list_path, 'w') as f:
            for file_path in self.alignment_files.split():
                f.write(f'{file_path}\n')
        
        try:
            concatenator = SequenceConcatenator(uppercase=False, printout=self.printout)
            num_taxa, total_length = concatenator.concatenate_files(
                file_list_path   = str(file_list_path),
                output_matrix    = self.s_matrix_dir,
                output_partition = self.s_model_dir
            )
            self.printout('metric', [('total_taxa', num_taxa), ('total_bp', f'{total_length}')])
            self.return_dict['supermatrix_taxa'] = num_taxa
            self.return_dict['supermatrix_bp']   = total_length
        except Exception as e:
            self.printout('error', f'Failed to create supermatrix: {str(e)}')
            sys.exit(1)
        finally:
            if file_list_path.exists():
                file_list_path.unlink()

    def iqtree(self):
        """
        Run IQ-TREE.
        IQ-TREE creates a tree from the supermatrix.
        """
        self.printout('metric', 'SuperMatrix')
        file_dir = (self.dir_super / 'SuperMatrix').resolve()
        s_matrix = (self.dir_super / 'super.matrix').resolve()
        s_model  = (self.dir_super / 'super.model').resolve()
        if self.super_bootstrap >= 1000:
            cmd = f'iqtree2 -s {s_matrix} -spp {s_model} -nt {self.threads} -bb {self.super_bootstrap} -m GTR+G -pre {file_dir} -redo'
        else:
            cmd = f'iqtree2 -s {s_matrix} -spp {s_model} -nt {self.threads} -m GTR+G -pre {file_dir} -redo'

        if self.subprocess_dir:
            result = run_logged_subprocess(cmd, self.subprocess_dir, 'IQ-TREE_SuperMatrix', shell=True, check=False)
        else:
            result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)

        if result.returncode != 0:
            self.printout('warning', 'IQ-TREE SuperMatrix failed (e.g. internal NNI assertion on large alignments). Retrying with fast tree only (-n 0).')
            cmd_fast = f'iqtree2 -s {s_matrix} -spp {s_model} -nt {self.threads} -n 0 -m GTR+G -pre {file_dir} -redo'
            if self.super_bootstrap >= 1000:
                cmd_fast = f'iqtree2 -s {s_matrix} -spp {s_model} -nt {self.threads} -n 0 -bb {self.super_bootstrap} -m GTR+G -pre {file_dir} -redo'
            if self.subprocess_dir:
                result = run_logged_subprocess(cmd_fast, self.subprocess_dir, 'IQ-TREE_SuperMatrix_fast', shell=True, check=False)
            else:
                result = subprocess.run(cmd_fast, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
            if result.returncode != 0:
                self.printout('error', 'IQ-TREE SuperMatrix failed')
                if result.stderr:
                    stderr_text = result.stderr.decode('utf-8') if isinstance(result.stderr, bytes) else str(result.stderr)
                    self.printout('error', stderr_text)
                sys.exit(1)

        if Path(self.supermatrix).exists():
            shutil.move(self.supermatrix, self.genetree)
            self._translate_supermatrix_tree_tips(self.genetree)

    def _translate_supermatrix_tree_tips(self, tree_path):
        if not self.renamed_map:
            return
        prefix_to_stem = {}
        for new_name, (_, source_file) in self.renamed_map.items():
            if '@' in new_name:
                prefix = new_name.split('@')[0]
            else:
                prefix = new_name.split('_')[0] if '_' in new_name else new_name
            stem = Path(source_file).stem
            if prefix not in prefix_to_stem:
                prefix_to_stem[prefix] = stem
        if not prefix_to_stem:
            return
        try:
            tree = ete3.Tree(tree_path)
            for node in tree.traverse():
                if node.is_leaf() and node.name:
                    node.name = prefix_to_stem.get(node.name, node.name)
            tree.write(format=1, outfile=tree_path)
        except Exception as e:
            self.printout('warning', f'Could not translate supermatrix tree tip labels: {e}')

    def bes(self):
        """
        Convert coalescent distance to molecular distance.
        """
        self.printout('metric', 'Convert Coalescent Distance to Molecular Distance')
        bes_runner        = BES(printout=self.printout)
        species_trees_dir = str(self.molecular_tree_path).replace('/SpeciesTree.molecular.tre', '')
        success           = bes_runner.run(self.finaltree, self.concat_tree, species_trees_dir)
        if not success:
            self.printout('error', 'BES failed to process trees')
            sys.exit(1)
        bes_trees_dir = Path(species_trees_dir) / 'BES_Trees'
        if bes_trees_dir.is_dir():
            self.return_dict['bes_tree_count'] = len(list(bes_trees_dir.iterdir()))