use std::env;
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand::rngs::StdRng;

extern crate regex;
use pest::Parser;

use newick_parser::node::{FlatTree, TraversalOrder};
use newick_parser::newick::{newick_to_tree, node_to_newick, NewickParser, Rule};

/// Removes a given leaf and its parent from the flat tree.
///
/// This function modifies the tree in place. It takes a flat tree and an index of a leaf to remove.
/// The corresponding leaf and its parent are removed, resulting in a correct phylogenetic tree with isolated nodes.
///
/// # Arguments
///
/// * `flat_tree` - A mutable reference to the flat tree.
/// * `index` - The index of the leaf to remove from the tree.
fn change_tree(flat_tree: &mut FlatTree, index: usize) {
    // Node indexes
    let parent_index = flat_tree[index]
        .parent
        .expect("The root is apparently a leaf.");
    let sister_index = if flat_tree[parent_index].left_child.unwrap() == index {
        flat_tree[parent_index].right_child.unwrap()
    } else {
        flat_tree[parent_index].left_child.unwrap()
    };

    let grandparent_index_opt = flat_tree[parent_index].parent;

    // The leaf and its parent are removed from the tree.
    flat_tree[parent_index].parent = None;

    // The sister of the leaf becomes the child of the grandparent.
    flat_tree[sister_index].parent = grandparent_index_opt;

    // Change the child of the grandparent from parent to sister.
    if let Some(grandparent_index) = grandparent_index_opt {
        if flat_tree[grandparent_index].left_child == Some(parent_index) {
            flat_tree[grandparent_index].left_child = Some(sister_index);
        } else {
            flat_tree[grandparent_index].right_child = Some(sister_index);
        }
    }
    // Depths and lengths will be recalculated later.
}

/// Removes all unsampled leaves and their parents from the flat tree.
///
/// This function takes a flat tree and a vector of indexes of leaves to remove. It removes all the corresponding leaves as well as their parents.
/// The function modifies the tree in place and does not return anything.
///
/// # Arguments
///
/// * `flat_tree` - A mutable reference to the flat tree.
/// * `list_indexes` - A vector of indexes of leaves to remove from the tree.
fn remove_all_unsampled(flat_tree: &mut FlatTree, list_indexes: &Vec<usize>) {
    for index in list_indexes {
        change_tree(flat_tree, *index);
    }
}

/// Extracts all leaves from the flat tree.
///
/// Leaves are nodes without children.
///
/// # Arguments
///
/// * `flat_tree` - A reference to the flat tree.
///
/// # Returns
///
/// A vector containing the indexes of all leaf nodes in the tree.
fn find_all_leaves(flat_tree: &FlatTree) -> Vec<usize> {
    flat_tree
        .iter(TraversalOrder::PreOrder)
        .enumerate()
        .filter(|(_, node)| node.left_child.is_none() && node.right_child.is_none())
        .map(|(i, _)| i)
        .collect()
}

/// Extracts all leaves that are present in the sampled species tree.
///
/// # Arguments
///
/// * `flat_tree` - A reference to the original flat tree.
/// * `flat_sampled_tree` - A reference to the flat tree of the sampled species tree.
///
/// # Returns
///
/// A vector containing the indexes of all leaves in the original tree that are also in the sampled tree.
fn find_all_extant_leaves(flat_tree: &FlatTree, flat_sampled_tree: &FlatTree) -> Vec<usize> {
    let sampled_names: Vec<String> = flat_sampled_tree
        .iter(TraversalOrder::PreOrder)
        .filter(|node| node.left_child.is_none() && node.right_child.is_none())
        .map(|node| node.name.clone())
        .collect();

    flat_tree
        .iter(TraversalOrder::PreOrder)
        .enumerate()
        .filter(|(_, node)| {
            node.left_child.is_none()
                && node.right_child.is_none()
                && sampled_names.contains(&node.name)
        })
        .map(|(i, _)| i)
        .collect()
}

/// Randomly samples a specified number of leaves from the tree.
///
/// This function only samples leaves that are present in both the original tree and the sampled species tree.
///
/// # Arguments
///
/// * `flat_tree` - A reference to the original flat tree.
/// * `flat_sampled_tree` - A reference to the flat tree of the sampled species tree.
/// * `n_sampled_nodes` - The number of leaves to sample.
/// * `rng` - A mutable reference to a random number generator.
///
/// # Returns
///
/// A vector containing the indexes of the sampled leaves.
fn sample_random_leaves(
    flat_tree: &FlatTree,
    flat_sampled_tree: &FlatTree,
    n_sampled_nodes: usize,
    rng: &mut StdRng,
) -> Vec<usize> {
    let leaves = find_all_extant_leaves(flat_tree, flat_sampled_tree);
    let sampled_leaves: Vec<usize> = leaves
        .choose_multiple(rng, n_sampled_nodes)
        .cloned()
        .collect();
    sampled_leaves
}

/// Determines which leaves should be removed from the tree.
///
/// This function calculates the complement of the sampled leaves within the set of all leaves, effectively identifying the leaves to be removed.
///
/// # Arguments
///
/// * `leaves` - A vector containing the indexes of all leaves in the tree.
/// * `sampled_leaves` - A vector containing the indexes of the sampled leaves.
///
/// # Returns
///
/// A vector containing the indexes of the leaves that are not in the sampled list.
fn leaves_to_be_removed(leaves: &Vec<usize>, sampled_leaves: &Vec<usize>) -> Vec<usize> {
    leaves
        .iter()
        .filter(|leaf| !sampled_leaves.contains(leaf))
        .cloned()
        .collect()
}

/// Finds the index of the root of the tree starting from a given leaf.
///
/// The function traverses up the tree from the given leaf node by following parent pointers until it reaches the root.
///
/// # Arguments
///
/// * `flat_tree` - A reference to the flat tree.
/// * `true_leaf` - The index of a leaf node from which to start the traversal.
///
/// # Returns
///
/// The index of the root node in the flat tree.
fn find_root(flat_tree: &FlatTree, true_leaf: usize) -> usize {
    let mut current_node = true_leaf;
    let mut current_parent = flat_tree[current_node].parent;
    while let Some(parent) = current_parent {
        current_node = parent;
        current_parent = flat_tree[current_node].parent;
    }
    current_node
}

/// Finds the leaves in the gene tree that correspond to the sampled species.
///
/// # Arguments
///
/// * `flat_gene_tree` - A reference to the flat tree of the gene tree.
/// * `leaf_names` - A vector of leaf names to find in the gene tree.
///
/// # Returns
///
/// A vector containing the indexes of the leaves in the gene tree that correspond to the given names.
fn find_leaves_in_gene_tree(flat_gene_tree: &FlatTree, leaf_names: &Vec<String>) -> Vec<usize> {
    flat_gene_tree
        .iter(TraversalOrder::PreOrder)
        .enumerate()
        .filter(|(_, node)| {
            node.left_child.is_none()
                && node.right_child.is_none()
                && leaf_names.contains(&node.name)
        })
        .map(|(i, _)| i)
        .collect()
}

/// Samples a gene tree by removing unsampled leaves and writes the result to a file.
///
/// # Arguments
///
/// * `sampled_leaves_names` - A vector of names of the sampled leaves.
/// * `leaves_to_be_removed_names` - A vector of names of the leaves to be removed.
/// * `gene_tree_path` - The path to the gene tree file.
/// * `gene_index` - The index of the gene tree.
/// * `output_dir` - The output directory where the sampled gene tree will be saved.
///
/// # Returns
///
/// A `Result` containing the Newick string of the sampled gene tree, or an `io::Error` if an error occurs.
fn one_gene_sample_to_string(
    sampled_leaves_names: &Vec<String>,
    leaves_to_be_removed_names: &Vec<String>,
    gene_tree_path: &PathBuf,
    gene_index: u32,
    output_dir: &str,
) -> Result<String, io::Error> {
    // Open the gene tree and convert it to a flat tree.
    let gene_tree_str = fs::read_to_string(&gene_tree_path)
        .unwrap_or_else(|err| panic!("Error reading file '{}': {}", gene_tree_path.to_string_lossy(), err));

    let gene_tree_str = gene_tree_str.trim();

    let mut pairs = NewickParser::parse(Rule::newick, gene_tree_str)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let mut node_tree = newick_to_tree(pairs.next().unwrap());

    let mut gene_tree = node_tree.pop().unwrap();
    gene_tree.zero_root_length();
    gene_tree.assign_depths(0.0);

    let mut flat_tree = gene_tree.to_flat_tree();

    // Find the indexes of the sampled leaves in the gene tree.
    let sampled_leaves = find_leaves_in_gene_tree(&flat_tree, sampled_leaves_names);
    let leaves_to_be_removed = find_leaves_in_gene_tree(&flat_tree, leaves_to_be_removed_names);

    // Remove unsampled leaves from the gene tree.
    remove_all_unsampled(&mut flat_tree, &leaves_to_be_removed);

    // Ensure the output directory exists
    fs::create_dir_all(output_dir)?;

    // Find the root of the new tree.
    let root_of_gene_tree = find_root(&flat_tree, sampled_leaves[0]);
    flat_tree.root = root_of_gene_tree;

    // Convert the flat tree back to a Node tree.
    let mut reconstructed_tree = flat_tree.to_node();

    // Update lengths based on depths.
    let root_depth = reconstructed_tree
        .depth
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Root depth not found"))?;
    reconstructed_tree.depths_to_lengths(root_depth);

    // Convert the gene tree to a Newick string.
    let reconstructed_newick = node_to_newick(&reconstructed_tree) + ";";

    // Save the gene tree as a Newick string in the output directory
    let gene_filename = format!("sampled_gene_{}.nwk", gene_index);
    let gene_filename = Path::new(output_dir).join(gene_filename);
    let mut gene_file = File::create(gene_filename)?;
    gene_file.write_all(reconstructed_newick.as_bytes())?;

    Ok(reconstructed_newick)
}

/// Samples the species tree and returns the Newick string along with sampled and removed leaf names.
///
/// The function performs the following steps:
/// 1. Reads and parses the species tree and the sampled species tree from files.
/// 2. Converts them to flat trees.
/// 3. Samples random leaves from the species tree that are present in the sampled species tree.
/// 4. Removes unsampled leaves from the species tree.
/// 5. Reconstructs the tree and updates node lengths based on depths.
/// 6. Converts the reconstructed tree to a Newick string and saves it to a file.
///
/// # Arguments
///
/// * `species_tree_path` - The path to the species tree file in Newick format.
/// * `sampled_species_tree_path` - The path to the sampled species tree file in Newick format.
/// * `n_sampled_nodes` - The number of leaves to sample.
/// * `output_dir` - The output directory where the sampled tree will be saved.
/// * `rng` - A mutable reference to a random number generator.
///
/// # Returns
///
/// A `Result` containing:
/// - The Newick string of the sampled species tree.
/// - A vector of names of sampled leaves.
/// - A vector of names of unsampled (removed) leaves.
///
/// If an error occurs, an `io::Error` is returned.
fn species_tree_sample_to_string(
    species_tree_path: &str,
    sampled_species_tree_path: &str,
    n_sampled_nodes: usize,
    output_dir: &str,
    rng: &mut StdRng,
) -> Result<(String, Vec<String>, Vec<String>), io::Error> {
    // Ensure the output directory exists
    let output_path = Path::new(output_dir);
    if !output_path.exists() {
        fs::create_dir_all(output_path)?;
    }

    // Read the species tree and the sampled species tree
    let species_tree_str = fs::read_to_string(species_tree_path)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let sampled_species_tree_str = fs::read_to_string(sampled_species_tree_path)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let species_tree_str = species_tree_str.trim();
    let sampled_species_tree_str = sampled_species_tree_str.trim();

    // Parse the species trees
    let mut pairs = NewickParser::parse(Rule::newick, species_tree_str)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let mut node_tree = newick_to_tree(pairs.next().unwrap());
    let mut species_tree = node_tree.pop().unwrap();

    let mut pairs_sampled = NewickParser::parse(Rule::newick, sampled_species_tree_str)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let mut node_tree_sampled = newick_to_tree(pairs_sampled.next().unwrap());
    let sampled_species_tree = node_tree_sampled.pop().unwrap();

    // Assign depths
    species_tree.zero_root_length();
    species_tree.assign_depths(0.0);
    let mut flat_tree = species_tree.to_flat_tree();

    let mut flat_sampled_tree = sampled_species_tree.to_flat_tree();

    // Sample random leaves
    let sampled_leaves = sample_random_leaves(&flat_tree, &flat_sampled_tree, n_sampled_nodes, rng);

    // Remove unsampled leaves
    let leaves = find_all_leaves(&flat_tree);
    let leaves_to_be_removed = leaves_to_be_removed(&leaves, &sampled_leaves);
    remove_all_unsampled(&mut flat_tree, &leaves_to_be_removed);

    // Update the root of the flat tree.
    let root_of_species_tree = find_root(&flat_tree, sampled_leaves[0]);
    flat_tree.root = root_of_species_tree;

    // Convert the flat tree back to a Node tree.
    let mut reconstructed_tree = flat_tree.to_node();

    // Update lengths based on depths.
    let root_depth = reconstructed_tree
        .depth
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Root depth not found"))?;
    reconstructed_tree.depths_to_lengths(root_depth);

    // Convert the species tree to a Newick string.
    let reconstructed_newick = node_to_newick(&reconstructed_tree) + ";";

    // Save the species tree as a Newick string in the output directory
    let species_filename = Path::new(output_dir).join("sampled_species_tree.nwk");
    let mut species_file = File::create(species_filename)?;
    species_file.write_all(reconstructed_newick.as_bytes())?;

    // Return the Newick string and the lists of sampled and removed leaf names.
    let sampled_leaves_names: Vec<String> = sampled_leaves
        .iter()
        .map(|i| flat_tree[*i].name.clone())
        .collect();
    let leaves_to_be_removed_names: Vec<String> = leaves_to_be_removed
        .iter()
        .map(|i| flat_tree[*i].name.clone())
        .collect();

    Ok((reconstructed_newick, sampled_leaves_names, leaves_to_be_removed_names))
}

/// Samples all gene trees in the specified range by removing unsampled leaves.
///
/// # Arguments
///
/// * `sampled_leaves_names` - A vector of names of the sampled leaves.
/// * `leaves_to_be_removed_names` - A vector of names of the leaves to be removed.
/// * `start_index` - The starting index of the gene trees to sample.
/// * `end_index` - The ending index of the gene trees to sample.
/// * `gene_trees_path` - The path to the directory containing the gene trees.
/// * `output_dir` - The output directory where the sampled gene trees will be saved.
fn sample_all_gene_trees(
    sampled_leaves_names: &Vec<String>,
    leaves_to_be_removed_names: &Vec<String>,
    start_index: usize,
    end_index: usize,
    gene_trees_path: &str,
    output_dir: &str,
) {
    let gene_trees_path = Path::new(gene_trees_path);
    for i in start_index..end_index {
        let gene_tree_filename = format!("gene_{}.nwk", i);
        let gene_tree_path = gene_trees_path.join("genes").join(gene_tree_filename);
        let _ = one_gene_sample_to_string(
            sampled_leaves_names,
            leaves_to_be_removed_names,
            &gene_tree_path,
            i as u32,
            &output_dir,
        );
    }
}

fn main() {
    // Read the arguments
    let args: Vec<String> = env::args().collect();
    // We use two trees, the species tree, and the sampled species tree, with leaves taken from the sampled species tree,
    // and from which we will sample leaves to obtain sampled gene trees. This serves for example if we have a complete tree and an
    // extant tree.
    // Ensure the correct number of arguments are provided
    if args.len() != 9 {
        eprintln!("Usage: {} <species_tree_path> <sampled_species_tree_path> <gene_trees_path> <n_sampled_nodes> <start_index> <end_index> <output_dir> <rng_seed>", args[0]);
        eprintln!("Received arguments: {:?}", args);
        return; // Exit early if the wrong number of arguments
    }

    let species_tree_path = &args[1];
    let sampled_species_tree_path = &args[2];
    let gene_trees_path = &args[3];
    let n_sampled_nodes = match args[4].parse::<usize>() {
        Ok(num) => num,
        Err(_) => {
            eprintln!("Error: n_sampled_nodes must be an integer. Received: {}", args[4]);
            eprintln!("All arguments: {:?}", args);
            return;
        }
    };
    let start_index = match args[5].parse::<usize>() {
        Ok(num) => num,
        Err(_) => {
            eprintln!("Error: start_index must be an integer. Received: {}", args[5]);
            eprintln!("All arguments: {:?}", args);
            return;
        }
    };
    let end_index = match args[6].parse::<usize>() {
        Ok(num) => num,
        Err(_) => {
            eprintln!("Error: end_index must be an integer. Received: {}", args[6]);
            eprintln!("All arguments: {:?}", args);
            return;
        }
    };
    let output_dir = &args[7];
    let rng_seed_str = &args[8];
    // Convert rng_seed_str to u64
    let seed = rng_seed_str.parse::<u64>().unwrap_or_else(|e| {
        eprintln!("Error parsing RNG seed: {}", e);
        std::process::exit(1);
    });

    let mut rng = StdRng::seed_from_u64(seed);
    // Sample the species tree
    let result = species_tree_sample_to_string(
        species_tree_path,
        sampled_species_tree_path,
        n_sampled_nodes,
        output_dir,
        &mut rng,
    );
    let (sampled_leaves_names, leaves_to_be_removed_names) = match result {
        Ok((_, sampled_names, removed_names)) => (sampled_names, removed_names),
        Err(e) => {
            eprintln!("Error during species tree sampling: {}", e);
            eprintln!("Species Tree Path: {}", species_tree_path);
            eprintln!("Number of Sampled Nodes: {}", n_sampled_nodes);
            eprintln!("Output Directory: {}", output_dir);
            return;
        }
    };

    // Sample the gene trees
    sample_all_gene_trees(
        &sampled_leaves_names,
        &leaves_to_be_removed_names,
        start_index,
        end_index,
        gene_trees_path,
        output_dir,
    );
}
