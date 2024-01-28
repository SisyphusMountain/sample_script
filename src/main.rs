/*
SCRIPT sample_rust
--------------------------------------------------------------------------------
This script takes:
<species_tree_path> <sampled_species_tree_path> <gene_trees_path> <n_sampled_nodes> <start_index> <end_index> <output_dir> <rng_seed>
- a tree in a format like such: ((a:1.0,b:1.0)c:1.0,d:2.0)e:0.0; (with necessary ; at the end, and no [&R] at the beginning)
- a path to the folder containing the corresponding gene trees in newick format, with filenames gene_folder/gene_i.nwk
- the number of leaves we want to keep on the species tree, and hence on all the gene trees.
- the start and index of the gene trees to which we want to apply the sampling. If you want to sample all
n gene trees you should put 0 n
- the output directory of the newick files
- the seed for the random number generators
--------------------------------------------------------------------------------
Output:
new files created: 
- output_folder/sampled_gene_i.nwk (sampled gene tree corresponding to the newick file gene_tree_path/gene_i.nwk)
- output_folder/sampled_species_tree.nwk (sampled species tree)
--------------------------------------------------------------------------------


This code should read a species tree and a list of gene trees, and sample a 
subset of the leaves.
To do this, we can do the same operation as in the horizontal gene transfer
simulation script, except that the new parent to the nodes removed now becomes
empty (because they are removed from the tree).
At the end, to find the new root of the tree, we take a leaf which was not
removed, and go upwards in the tree to find the new root.
--------------------------------
STEPS:
1. Open the species tree.
(With newick_to_tree)
2. Convert it to a flat_tree object
(With node_to_flat)
3. Sample a subset of the leaves, and save it to a vector.
(With sample_random_leaves)
4. Reconstruct the sampled species tree, by SPR moves, disconnecting unsampled leaves
(Use change_tree to remove any given node from the flat_tree)
    4.1 Find the new root of the tree by taking an unsampled leaf and going upwards in the tree.
    (With find_root)

5. Save the newick string of the sampled species tree to a file.
6. For each gene newick file:
    6.1. Open the gene tree.
    6.2. Convert it to a flat_tree object.
        6.2.1. We have to find the nodes in the gene tree which correspond to the nodes in the species tree.
        Since the leaves should all have names, we can use the names to find the corresponding nodes.
        We just need to extract all names of sampled nodes, and then find the corresponding node indexes in the gene tree.
        (With find_leaves_in_gene_tree)
    6.3. Reconstruct the sampled gene tree, by SPR moves, disconnecting unsampled leaves, using the leaves_in_gene_tree vector
    6.4. Save the newick string of the sampled gene tree to a file.



The script samples leaves, but only from the extant ones.
*/

#[macro_use]
extern crate pest_derive;
#[allow(dead_code)]
use pest::Parser;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::env;
use std::fs::File;
use rand::seq::SliceRandom;
extern crate regex;
use rand::SeedableRng;
use rand::rngs::StdRng;
#[derive(Clone, Debug)]


struct Node {
    name: String,
    left_child: Option<Box<Node>>,
    right_child: Option<Box<Node>>,
    parent: Option<usize>, // Using index for parent
    depth: Option<f64>,
    length: f64,
}
#[derive(Clone, Debug)]
struct FlatNode {
    name: String,
    left_child: Option<usize>,
    right_child: Option<usize>,
    parent: Option<usize>,
    depth: Option<f64>,
    length: f64,
}
#[derive(Parser)]
#[grammar = "newick.pest"]
struct NewickParser;

// --------------------------------

fn give_depth(node: &mut Node, depth: f64) {
    /*
    Recursively gives the depth of nodes in a "Node" arborescent tree.
    ----------------------------------
    Input:
    - A node, and its depth in the tree
    ----------------------------------
    Output:
    - None. The tree is modified in place via a mutable reference (&mut)
    */
    node.depth = Some(depth);
    if let Some(left_child) = &mut node.left_child {
        give_depth(left_child, depth + left_child.length);
    }
    if let Some(right_child) = &mut node.right_child {
        give_depth(right_child, depth + right_child.length);
    }
}
// The following two functions are used to convert a newick string into an arborescent Tree structure, where each "Node" object owns its two children in a Box object.
fn newick_to_tree(pair: pest::iterators::Pair<Rule>) -> Vec<Node> {
    /*
    The input .nwk file of the script is a concatenation of newick trees, separated by ";".
    This function returns a vector containing each of these trees.
    --------------------------------
    Input:
    - A "pair" object of type "newick" (see newick.pest grammar) which is a newick string representing one single tree.
    ----------------------------------
    Output:

    */
    let mut vec_trees:Vec<Node> = Vec::new();
    for inner in pair.into_inner() {
        let tree = handle_pair(inner);
        vec_trees.push(tree.unwrap());
    }
    return vec_trees
}
fn handle_pair(pair: pest::iterators::Pair<Rule>) -> Option<Node> {
    /* 
    Recursively parses a newick string representing a single newick tree, according to the grammar contained in newick.pest.
    ----------------------------------
    Input:
    - A "pair" object representing a part of this string.
    ----------------------------------
    Output:
    - Either a node if the input string represents a node, or nothing in cases where it does not (e.g. if the string is the length of a node).
    */
    match pair.as_rule() {
        Rule::newick => {
            // newick = { subtree ~ ";" }
            // For now, we only do the implementation for one single tree.
            // This case is trivial : we just call handle_pair on the only tree we found.
            for inner in pair.into_inner() {
                if inner.as_rule() == Rule::subtree {
                    return handle_pair(inner);
                }
            }
            None
        },
        Rule::subtree => {
            // subtree = { leaf | internal }
            // Subtree is like the choice between an inner node and a leaf. Either way, we need to pass it to handle_pair.
            handle_pair(pair.into_inner().next().unwrap())
        },
        Rule::leaf => {
            // Choose default values for the name and length of the leaf.
            // The defaults are an empty string and a 0.0 length.
            let mut name = String::new();
            let mut length = 0.0;

            // leaf = { NAME? ~ ":"? ~ LENGTH? }
            // Therefore, we use a match statement to handle the cases NAME and LENGTH because ":" is unimportant.
            for inner_pair in pair.into_inner() {
                match inner_pair.as_rule() {
                    Rule::NAME => {
                        name = inner_pair.as_str().to_string();
                    }
                    Rule::LENGTH => {
                        let val = inner_pair.as_str();
                        if let Err(_) = val.parse::<f64>() {
                            println!("Failed to parse LENGTH: {}", val);
                        }
                        length = val.parse::<f64>().unwrap_or(0.0);
                    },
                    
                    _ => {} // Ignore other rules
                }
            }
            let node = Node {
                name: name,
                left_child: None,
                right_child: None,
                parent: None,
                depth: None,
                length: length,
            };
            Some(node)
        },
        Rule::internal => {
            // internal = { "(" ~ subtree ~ "," ~ subtree ~ ")" ~ NAME? ~ ":"? ~ LENGTH? }
            
            // Initialize default values for the name and length of the internal node.
            // The defaults are an empty string and a 0.0 length.
            let mut name = String::new();
            let mut length = 0.0;
        
            let mut first_subtree = None;
            let mut second_subtree = None;
        
            // Iterate through the inner rules without assuming their order
            for inner_pair in pair.into_inner() {
                match inner_pair.as_rule() {
                    Rule::subtree => {
                        let subtree = handle_pair(inner_pair).unwrap();
                        if first_subtree.is_none() {
                            first_subtree = Some(subtree);
                        } else {
                            second_subtree = Some(subtree);
                        }
                    },
                    Rule::NAME => {
                        name = inner_pair.as_str().to_string();
                    },
                    Rule::LENGTH => {
                        let val = inner_pair.as_str();
                        if let Err(_) = val.parse::<f64>() {
                            println!("Failed to parse LENGTH: {}", val);
                        }
                        length = val.parse::<f64>().unwrap_or(0.0);
                    },
                    
                    _ => {} // Ignore other rules
                }
            }
        
            let node = Node {
                name,
                left_child: first_subtree.map(Box::new),
                right_child: second_subtree.map(Box::new),
                parent: None,
                depth: None,
                length,
            };
            Some(node)
        }
        Rule::NAME | Rule::LENGTH => {
            // We should never directly handle these outside their containing rules.
            None
        },
    }
}

// Convert a Tree in "Node" form to a Newick string.
fn node_to_newick(node: &Node) -> String {
    /* Takes a node and returns the corresponding subtree in Newick format.
        --------------------------------
        INPUT:
            - node: the node to convert to Newick format.
        OUTPUT:
            - the Newick representation of the subtree rooted at node.
        Warning: rounds the lengths to 6 decimal places.
    */
    if let (Some(left_child), Some(right_child)) = (&node.left_child, &node.right_child) {
        // This is an internal node with both left and right children.
        format!(
            "({},{}){}:{:.6}",
            node_to_newick(left_child),
            node_to_newick(right_child),
            node.name,
            node.length
        )
    } else {
        // This is a leaf node.
        format!("{}:{:.6}", node.name, node.length)
    }
}
// Convert from FlatNode to Node
fn flat_to_node(flat_tree: &[FlatNode], index: usize, parent_index: Option<usize>) -> Option<Node> {
    /*
    This function converts a flat tree into an arborescent tree recursively.
    To use it, give the flat_tree, as well as the index of the root. The vector will be traversed recursively, 
    following the descendants of each node being examined.
    ----------------------------------
    Input: 
    - A flat tree, the index of a node in the flat tree, and the index of its parent.
    ----------------------------------
    Output:
    - The corresponding node in the arborescent tree.
    ----------------------------------
    Warning: can bug if applied directly to a non-root node, which will be mistaken for a root node.
    */
    let flat_node = &flat_tree[index];
    let left_child = flat_node.left_child.and_then(|i| {
        flat_to_node(flat_tree, i, Some(index)).map(Box::new)
    });
    let right_child = flat_node.right_child.and_then(|i| {
        flat_to_node(flat_tree, i, Some(index)).map(Box::new)});
    
    Some(Node {
        name: flat_node.name.clone(),
        left_child: left_child,
        right_child: right_child,
        parent: parent_index,
        depth: flat_node.depth,
        length: flat_node.length,
    })
}
// Convert from Node to FlatNode
fn node_to_flat(node: &Node, flat_tree: &mut Vec<FlatNode>, parent: Option<usize>) -> usize {
    /* Transforms the arborescent tree into a "linear tree", which is just a vector of nodes with parent and descendants.
    ----------------------------------
    Input:
    - The root node, which contains the whole arborescent tree
    - The flat_tree vector, which is to be filled by the function
    ----------------------------------
    Output:
    - The index of the node currently being added to flat_tree (usize because the indexing of Rust vectors can only be usize).
    */
    let index = flat_tree.len();
    flat_tree.push(FlatNode {
        name: node.name.clone(),
        left_child: None,  // Will fill this in a moment
        right_child: None, // Will fill this in a moment
        parent: parent,
        depth: node.depth,
        length: node.length,
    });

    if let Some(left) = &node.left_child {
        let left_index = node_to_flat(left, flat_tree, Some(index));
        flat_tree[index].left_child = Some(left_index);
    }

    if let Some(right) = &node.right_child {
        let right_index = node_to_flat(right, flat_tree, Some(index));
        flat_tree[index].right_child = Some(right_index);
    }

    index
}
// Remove any given leaf from the tree
fn change_tree(flat_tree: &mut Vec<FlatNode>, index: usize){
    /* 
    Takes a flat tree and an index of leaf to remove from the flat tree.
    Removes the corresponding leaf as well as its parent.
    Each application of this function should take a correct phylogenetic tree,
    and output an object representing a correct phylogenetic tree along with isolated nodes.
    ----------------------------------
    Input:
    - A flat tree, and the index of a leaf to remove from the tree.
    ----------------------------------
    Output:
    - The flat tree, with the leaf removed.
    */

    // --------------------------------
    // Node indexes
    let parent_index = flat_tree[index].parent.expect("The root is apparently a leaf.");
    let sister_index;
    if flat_tree[parent_index].left_child.unwrap() == index {
        sister_index = flat_tree[parent_index].right_child.unwrap()} else 
        {sister_index = flat_tree[parent_index].left_child.unwrap()};

    let grandparent_index_opt = flat_tree[parent_index].parent;




        // The leave and its parent are removed from the tree.
        flat_tree[parent_index].parent = None;

        // The sister of the leaf becomes the child of the grandparent, because its parent has lost one of its child and hence disappears.
        flat_tree[sister_index].parent = grandparent_index_opt;


        // Changing attributes of the parent of the removed leaf is unimportant, since it will never be referenced again, and hence
        // will not be seen in a tree traversal from the root.

        // Changing the child of the grandparent which is the parent to the sister.
        if let Some(grandparent_index) = grandparent_index_opt {
            if flat_tree[grandparent_index].left_child == Some(parent_index) {
                flat_tree[grandparent_index].left_child = Some(sister_index)} else
                {flat_tree[grandparent_index].right_child = Some(sister_index)
                };
            }
        // None of the depths have to be updated. The lengths are not updated either, since they will
        // be reconstructed from the depths at the end.
    }
// Apply change_tree to remove all unsampled leaves
fn remove_all_unsampled(flat_tree: &mut Vec<FlatNode>, list_indexes: &Vec<usize>) {
    /*
    Takes a flat tree and a vector of indexes of leaves to remove from the tree.
    Removes all the corresponding leaves as well as their parents.
    Acts in place and does not return anything.
    ----------------------------------
    Input:
    - A flat tree, and a vector of indexes of leaves to remove from the tree.
    ----------------------------------
    Output:
    - The flat tree, with the leaves removed.
    */
    for index in list_indexes {
        change_tree(flat_tree, *index);
    }
}
// Extract leaves from the tree
fn find_all_leaves(flat_tree: &Vec<FlatNode>) -> Vec<usize> {
    // Leaves are exactly nodes without children, so they are easy to find.
    let mut leaves: Vec<usize> = Vec::new();
    for (i, node) in flat_tree.iter().enumerate() {
        if node.left_child.is_none() && node.right_child.is_none() {
            leaves.push(i);
        }
    }
    leaves
}

fn find_all_extant_leaves(flat_tree: &Vec<FlatNode>, flat_sampled_tree: &Vec<FlatNode>) -> Vec<usize> {
    let mut leaves: Vec<usize> = Vec::new();
    for (i, node) in flat_tree.iter().enumerate() {
        if node.left_child.is_none() && node.right_child.is_none() {
            // Check if this leaf node is present in flat_sampled_tree
            if flat_sampled_tree.iter().any(|sampled_node| sampled_node.name == node.name) {
                leaves.push(i);
            }
        }
    }
    leaves
}
// Choose random leaves uniformly
fn sample_random_leaves(flat_tree: &Vec<FlatNode>, flat_sampled_tree: &Vec<FlatNode>, n_sampled_nodes: usize, rng: &mut StdRng) -> Vec<usize> {
    // This function samples n_sampled_nodes random leaves from the tree.
    // We make it so that it only samples extant leaves
    // First, we find all the leaves.
    let leaves = find_all_extant_leaves(flat_tree, flat_sampled_tree);
    // Then, we sample n_sampled_nodes of them.
    let sampled_leaves: Vec<usize> = leaves.choose_multiple(rng, n_sampled_nodes).cloned().collect();
    sampled_leaves
}
// Find the leaves to be removed (Complement of leaves to be kept)
fn leaves_to_be_removed(leaves: &Vec<usize>, sampled_leaves: &Vec<usize>) -> Vec<usize> {
    // This function takes a vector of all the leaves in the tree, and a vector of the sampled leaves.
    // It returns a vector of the indexes of the leaves to be removed.
    let mut leaves_to_be_removed: Vec<usize> = Vec::new();
    for leaf in leaves {
        if !sampled_leaves.contains(leaf) {
            leaves_to_be_removed.push(*leaf);
        }
    }
    leaves_to_be_removed
}
// Give the index of the root of the tree
fn find_root(flat_tree: &Vec<FlatNode>, true_leaf: usize) -> usize {
    /*
    Uses an index of a leaf to find the root of the tree.
     */
    let mut current_node = true_leaf;
    let mut current_parent = flat_tree[current_node].parent;
    while current_parent.is_some() {
        current_node = current_parent.unwrap();
        current_parent = flat_tree[current_node].parent;
    }
    current_node
}
// Find the leaves in the gene tree corresponding to the leaves in the species tree
fn find_leaves_in_gene_tree(flat_gene_tree: &Vec<FlatNode>, leave_names: &Vec<String>) -> Vec<usize> {
    /*
    Takes a gene tree and a vector of names of leaves in the species tree.
    Returns a vector of indexes of the leaves in the gene tree.
    ----------------------------------
    Input:
    - A gene tree, and a vector of names of leaves in the species tree.
    ----------------------------------
    Output:
    - A vector of indexes of the leaves in the gene tree.
    */
    let mut leaves_in_gene_tree: Vec<usize> = Vec::new();
    for (i, node) in flat_gene_tree.iter().enumerate() {
        if node.left_child.is_none() && node.right_child.is_none() {
            if leave_names.contains(&node.name) {
                leaves_in_gene_tree.push(i);
            }
        }
    }
    leaves_in_gene_tree
}
// Corrects lengths of the tree after removing some nodes
fn depths_to_lengths(node: &mut Node, parent_depth: f64) {
    /*
    Computes the lengths of nodes from the depths.
    Necessary because the gene transfer functions only modify the depth of the
    transfered node and not lengths for convenience.
    ----------------------------------
    Input:
    - A "Node" object representing an arborescent tree.
    - The depth of its parent
    ----------------------------------
    Output:
    - None. The tree is modified in place via a mutable reference.
    */
    let depth = node.depth.unwrap();
    if node.parent.is_some() {
        node.length = depth - parent_depth;
    }
    if let Some(left_child) = &mut node.left_child {
        depths_to_lengths(left_child, depth);
    }
    if let Some(right_child) = &mut node.right_child {
        depths_to_lengths(right_child, depth);
    }
}
// Returns and saves a gene tree with sampled leaves
fn one_gene_sample_to_string(sampled_leaves_names: &Vec<String>, leaves_to_be_removed_names: &Vec<String>, gene_tree_path: &PathBuf, gene_index: u32, output_dir: &str) -> Result<String, io::Error> {
    /*
    Samples all the species in the gene tree.
    1. Find the indexes of the sampled leaves thanks to their names.
    2. Construct the gene tree by removing the unsampled leaves from the species tree.
        2.1. Apply the function remove_all_unsampled to the flat tree.
        2.2. Use any of the leaves in sampled_leaves to find the root of the new tree.
        2.3. Convert the flat tree to a classical tree by using the root index. This will effectively ignore all the nodes which are now unsampled.
    3. Convert the gene tree to a Newick string.
    ----------------------------------
    Input:
    - The path to the gene tree file.
    - The number of sampled nodes.
    - The index of the gene tree.
    - The path to the output directory.
    ----------------------------------
    Output:
    - The Newick string of the sampled gene tree.
        */
    // 0. Open the gene tree and convert it to a flat tree.
    //println!("gene_tree_path {:?}", gene_tree_path);
    let gene_tree = fs::read_to_string(&gene_tree_path)
        .unwrap_or_else(|err| panic!("Error reading file '{}': {}", gene_tree_path.to_string_lossy(), err));


    let gene_tree = gene_tree.trim();
    let mut gene_tree = NewickParser::parse(Rule::newick, gene_tree)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let mut gene_tree = newick_to_tree(gene_tree.next().unwrap());
    let root_length = gene_tree[0].length.clone();
    let _ = give_depth(&mut gene_tree[0], root_length);
    let mut copied_flat_tree: Vec<FlatNode> = Vec::new();
    // Populate the flat tree with the gene tree.
    let _ = node_to_flat(&gene_tree[0], &mut copied_flat_tree, None);
    // 1. Find the indexes of the removed and sampled leaves thanks to their names.
    let sampled_leaves = find_leaves_in_gene_tree(&copied_flat_tree, sampled_leaves_names);
    let leaves_to_be_removed = find_leaves_in_gene_tree(&copied_flat_tree, leaves_to_be_removed_names);
    // 2. Construct the gene tree by removing the unsampled leaves from the species tree.
    // 2.1. Apply the function remove_all_unsampled to the flat tree.
    remove_all_unsampled(&mut copied_flat_tree, &leaves_to_be_removed);
    // Ensure the output directory exists
    fs::create_dir_all(output_dir)?;
    // 2.2. Use any of the leaves in sampled_leaves to find the root of the new tree.
    let root_of_gene_tree = find_root(&copied_flat_tree, sampled_leaves[0]);
    // 2.3. Convert the flat tree to a classical tree by using the root index. This will effectively ignore all the nodes which are now unsampled.
    let mut reconstructed_tree = flat_to_node(&copied_flat_tree, root_of_gene_tree, None)
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Failed to convert flat tree to node tree"))?;

    // The lengths of the nodes are not correct, so we need to update them.
    let root_depth = reconstructed_tree.depth
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Root depth not found"))?;
    depths_to_lengths(&mut reconstructed_tree, root_depth);

    // 3. Convert the gene tree to a Newick string.
    let reconstructed_newick = node_to_newick(&reconstructed_tree) + ";";

    // Save the gene tree as a Newick string in the genes directory
    let gene_filename = format!("sampled_gene_{}.nwk", gene_index);
    let gene_filename = Path::new(output_dir).join(gene_filename);
    // println!("Writing to file: {:?}", gene_filename);
    let mut gene_file = File::create(gene_filename)?;
    gene_file.write_all(reconstructed_newick.as_bytes())?;

    Ok(reconstructed_newick)
}
fn species_tree_sample_to_string(species_tree_path: &str, sampled_species_tree_path: &str, n_sampled_nodes: usize, output_dir: &str, rng: &mut StdRng) -> Result<(String, Vec<String>, Vec<String>), io::Error> {
    /*
    Samples the given number of nodes in the species tree, and returns the list of names of sampled nodes, as well as unsampled nodes.
    1. Convert the species tree to a flat tree.
    2. Sample the leaves.
    3. Reconstruct the sampled species tree, by SPR moves, disconnecting unsampled leaves.
        3.1. Apply the function remove_all_unsampled to the flat tree.
        3.2. Use any of the leaves in sampled_leaves to find the root of the new tree.
        3.3. Convert the flat tree to a classical tree by using the root index. This will effectively ignore all the nodes which are now unsampled.
    4. Convert the species tree to a Newick string.
    ----------------------------------
    Input:
    - The path to the species tree file.
    - The number of sampled nodes.
    - The path to the output directory.
    ----------------------------------
    Output:
    - The Newick string of the sampled species tree.
    - The list of names of sampled nodes.
    - The list of names of unsampled nodes.
    */
    // -1. Ensure the output directory exists
    let output_path = Path::new(output_dir);
    if !output_path.exists() {
        fs::create_dir_all(output_path)?;
    }

    // 0. Open the species tree and convert it to a flat tree.
    let species_tree = fs::read_to_string(species_tree_path)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let sampled_species_tree = fs::read_to_string(sampled_species_tree_path)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let species_tree = species_tree.trim();
    let sampled_species_tree = sampled_species_tree.trim();
    let mut species_tree = NewickParser::parse(Rule::newick, species_tree)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let mut sampled_species_tree = NewickParser::parse(Rule::newick, sampled_species_tree)
        .map_err(|e| io::Error::new(io::ErrorKind::Other, e.to_string()))?;
    let mut species_tree = newick_to_tree(species_tree.next().unwrap());
    let sampled_species_tree = newick_to_tree(sampled_species_tree.next().unwrap());
    let root_length = species_tree[0].length.clone();
    let _ = give_depth(&mut species_tree[0], root_length);
    let mut flat_tree: Vec<FlatNode> = Vec::new();
    // Populate the flat tree with the species tree.
    let _ = node_to_flat(&species_tree[0], &mut flat_tree, None);
    let mut flat_sampled_tree = Vec::new();
    let _ = node_to_flat(&sampled_species_tree[0], &mut flat_sampled_tree, None);

    // 1. Sample the leaves.
    let sampled_leaves = sample_random_leaves(&flat_tree, &flat_sampled_tree, n_sampled_nodes, rng);
    // 2. Construct the species tree by removing the unsampled leaves from the species tree.
    // 2.1. Apply the function remove_all_unsampled to the flat tree.
    let leaves_to_be_removed = leaves_to_be_removed(&find_all_leaves(&flat_tree), &sampled_leaves);
    remove_all_unsampled(&mut flat_tree, &leaves_to_be_removed);
    // 2.2. Use any of the leaves in sampled_leaves to find the root of the new tree.
    let root_of_species_tree = find_root(&flat_tree, sampled_leaves[0]);
    // 2.3. Convert the flat tree to a classical tree by using the root index. This will effectively ignore all the nodes which are now unsampled.
    let mut reconstructed_tree = flat_to_node(&flat_tree, root_of_species_tree, None)
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Failed to convert flat tree to node tree"))?;
    
    // The lengths of the nodes are not correct, so we need to update them.
    let root_depth = reconstructed_tree.depth
        .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "Root depth not found"))?;
    depths_to_lengths(&mut reconstructed_tree, root_depth);

    // 3. Convert the species tree to a Newick string.
    let reconstructed_newick = node_to_newick(&reconstructed_tree) + ";";

    // Save the species tree as a Newick string in the output directory
    let species_filename = format!("sampled_species_tree.nwk");
    let species_filename = Path::new(output_dir).join(species_filename);
    let mut species_file = File::create(species_filename)?;
    species_file.write_all(reconstructed_newick.as_bytes())?;

    // 4. Return the Newick string of the sampled species tree, as well as the list of names of sampled nodes, and the list of names of unsampled nodes.
    let sampled_leaves_names: Vec<String> = sampled_leaves.iter().map(|i| flat_tree[*i].name.clone()).collect();
    let leaves_to_be_removed_names: Vec<String> = leaves_to_be_removed.iter().map(|i| flat_tree[*i].name.clone()).collect();
    Ok((reconstructed_newick, sampled_leaves_names, leaves_to_be_removed_names))
}
fn sample_all_gene_trees(sampled_leaves_names: &Vec<String>, leaves_to_be_removed_names: &Vec<String>, start_index: usize, end_index: usize, gene_trees_path: &str, output_dir: &str) {
    /*
    For each of the gene trees in the specified range, open the corresponding files, and sample the leaves.
    ----------------------------------
    Input:
    - The number of sampled nodes.
    - The index of the first gene tree to sample.
    - The index of the last gene tree to sample.
     */
    let gene_trees_path = Path::new(gene_trees_path);
    for i in start_index..end_index {
        let gene_tree_filename = format!("gene_{}.nwk", i);
        let gene_tree_path = gene_trees_path.join("genes").join(gene_tree_filename);
        let _ = one_gene_sample_to_string(
            sampled_leaves_names,
            leaves_to_be_removed_names,
            &gene_tree_path,
            i as u32,
            &output_dir
        );
    }
}

fn main() {
    // Read the arguments
    let args: Vec<String> = env::args().collect();

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
    let result = species_tree_sample_to_string(species_tree_path, sampled_species_tree_path, n_sampled_nodes, output_dir, &mut rng);
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


    sample_all_gene_trees(&sampled_leaves_names, &leaves_to_be_removed_names, start_index, end_index, gene_trees_path, output_dir);
}
