# Introduction
This is a Rust program designed to sample a certain subset of leaves uniformly in a species tree as well as all the corresponding gene trees.




# Installation
Building a project written in Rust is very easy.
1. Install the [Cargo package manager](https://doc.rust-lang.org/cargo/getting-started/installation.html).
2. Clone the repository.
3. Navigate inside the directory in the terminal.
4. Run the command ```cargo build --release```. This will automatically download all necessary dependencies and compile the binary inside a new target directory.

# Algorithm
This Rust script chooses species uniformly among leaves of the species tree, then removes all other leaves. It then removes the leaves with the same names from the gene trees chosen.


# Usage
Input:
1. A .nwk file representing the species tree in a format like 
```((a:1.0,b:1.0)c:1.0,d:2.0)e:0.0;```.
2. A path for the folder ```gene_folder``` containing gene trees, which should be named ```gene_folder/gene_i.nwk```.
3. An integer giving the number of leaves to keep on the sampled trees.
4. The start index ```i``` of the first gene to sample (inclusive).
5. The end index ```j``` of the last gene to sample (exclusive).
6. The output directory ```output_dir```.

```Usage: target/release/sample_script <species_tree_path> <gene_trees_path> <n_sampled_nodes> <start_index> <end_index> <output_dir> <seed>```

Output:
```output_dir``` will contain the sampled species tree ```output_dir/sampled_species_tree.nwk``` as well as the sampled gene trees ```output_dir/sampled_gene_i.nwk```.
