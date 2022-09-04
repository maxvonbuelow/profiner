#pragma once
#include <vector>
#include <thread>
#include <iostream>
#include <bitset>
#include <queue>
#include <array>

/*
 * Helper library for computing bytes:flops ratios
 * (binning reuse distance)
 *
 * By Scott Pakin <pakin@lanl.gov>
 *    Rob Aulwes <rta@lanl.gov>
 */
// https://github.com/lanl/Byfl/blob/master/LICENSE.md

// #include "byfl.h"
#include <unordered_map>

// An RDnode is one node in a reuse-distance tree.
class RDnode {
private:
  RDnode* left;         // Left child
  RDnode* right;        // Right child

  // Fix the node's weight (subtree size).
  inline void fix_node_weight();

  // Fix the weight of all nodes along the path to a given time.
  inline void fix_path_weights(uint64_t time);

  // Splay a value to the top of the tree, returning the new tree.
  inline RDnode* splay(uint64_t target);

public:
  uint64_t address;     // Address from trace
  uint64_t time;        // Time of the address's last access
  uint64_t weight;      // Number of items in this subtree (self included)

  // Initialize a new RDnode with a given address and timestamp (defaulting to
  // dummy values).
  inline RDnode();
  inline RDnode(uint64_t address, uint64_t time);

  // Reinitialize an existing RDnode with a given address and timestamp.
  inline void initialize(uint64_t address, uint64_t time);

  // Insert a node into the tree and return the new tree.
  inline RDnode* insert(RDnode* new_node);

  // Remove a timestamp from the tree and return the new tree and the node that
  // was deleted.
  inline RDnode* remove(uint64_t timestamp, RDnode** removed_node);

  // Remove all timestamps less than a given value from the tree and from a
  // given histogram, and return the new tree.
//   RDnode* prune_tree(uint64_t timestamp, addr_to_time_t* histogram);

  // Return the number of nodes in a splay tree whose timestamp is larger than
  // a given value.
  inline uint64_t tree_dist(uint64_t timestamp);
};


// fix_node_weight() sets the weight of a given node to the sum of its
// immediate children's weight plus one.
void RDnode::fix_node_weight()
{
  uint64_t new_weight = 1;
  if (left != nullptr)
    new_weight += left->weight;
  if (right != nullptr)
    new_weight += right->weight;
  weight = new_weight;
}


// fix_path_weights() fixes node weights along the path to a given time.
void RDnode::fix_path_weights(uint64_t target)
{
  // Do an ordinary binary tree search for target -- which we expect not to
  // find -- but change child pointers to parent pointers as we go (instead of
  // requiring extra memory to maintain our path back to the root).
  RDnode* parent = nullptr;
  RDnode* node = this;
  while (node != nullptr) {
    RDnode* child;
    if (target < node->time) {
      child = node->left;
      node->left = parent;
    }
    else {
      child = node->right;
      node->right = parent;
    }
    parent = node;
    node = child;
  }

  // Walk back up the tree, fixing weights and child pointers as we go.
  while (parent != nullptr) {
    RDnode* prev_node = node;
    node = parent;
    if (target < node->time) {
      // We borrowed our left child's pointer.
      parent = node->left;
      node->left = prev_node;
    }
    else {
      // We borrowed our right child's pointer.
      parent = node->right;
      node->right = prev_node;
    }
    node->fix_node_weight();
  }
}


// splay() splays a value (or a nearby value if the value doesn't appear in the
// tree) to the top of the tree, returning the new tree.
RDnode* RDnode::splay(uint64_t target)
{
  RDnode* node = this;
  RDnode new_node;
  new_node.left = nullptr;
  new_node.right = nullptr;
  RDnode* left = &new_node;
  RDnode* right = &new_node;

  while (true) {
    if (target < node->time) {
      if (node->left == nullptr)
        break;
      if (target < node->left->time) {
        // Rotate right
        RDnode* parent = node->left;
        node->left = parent->right;
        parent->right = node;
        node = parent;

        // Fix weights.
        node->right->fix_node_weight();
        node->fix_node_weight();
        if (node->left == nullptr)
          break;
      }

      // Link right
      right->left = node;
      right = node;
      node = node->left;
    }
    else
      if (target > node->time) {
        if (node->right == nullptr)
          break;
        if (target > node->right->time) {
          // Rotate left
          RDnode* parent = node->right;
          node->right = parent->left;
          parent->left = node;
          node = parent;

          // Fix weights.
          node->left->fix_node_weight();
          node->fix_node_weight();
          if (node->right == nullptr)
            break;
        }

        // Link left
        left->right = node;
        left = node;
        node = node->right;
      }
      else
        break;
  }

  // Assemble the final tree.
  left->right = node->left;
  right->left = node->right;
  node->left = new_node.right;
  node->right = new_node.left;

  // Fix weights up to the node from its previous position.
  if (node->left != nullptr)
    node->left->fix_path_weights(node->time);
  if (node->right != nullptr)
    node->right->fix_path_weights(node->time);
  return node;
}

// insert() inserts a new node into a splay tree and returns the new tree.
// Duplicates insertions produce undefined behavior.  Insertions into NULL
// trees produce undefined behavior.  (The caller should check for the
// first-insertion and allocate memory accordingly.)
RDnode* RDnode::insert(RDnode* new_node)
{
  // Handle some simple cases.
  RDnode* node = this;
  node = node->splay(new_node->time);
  if (new_node->time == node->time)
    // The timestamp is already in the tree.  This should never happen when the
    // tree is used for reuse-distance calculations.
    abort();

  // Handle the normal cases.
  if (new_node->time > node->time) {
    new_node->right = node->right;
    new_node->left = node;
    node->right = nullptr;
  }
  else {
    new_node->left = node->left;
    new_node->right = node;
    node->left = nullptr;
  }
  node->fix_node_weight();
  new_node->fix_node_weight();
  return new_node;
}


// remove() deletes a timestamp from the tree and returns the new tree and the
// deleted node.  Missing timestamps produce undefined behavior.
RDnode* RDnode::remove(uint64_t target, RDnode** removed_node)
{
  RDnode* node = this;
  node = node->splay(target);
  if (node->time != target)
    // Not found
    abort();
  RDnode* new_root;
  if (node->left == nullptr)
    // Smallest value in the tree
    new_root = node->right;
  else {
    // Any other value
    new_root = node->left->splay(target);
    if (new_root != nullptr) {
      new_root->right = node->right;
      if (new_root->right != nullptr)
        new_root->right->fix_node_weight();
      new_root->fix_node_weight();
    }
  }
  *removed_node = node;
  return new_root;
}


// Remove all timestamps less than a given value from the tree and from a given
// histogram, and return the new tree and new set of symbols.
// RDnode* RDnode::prune_tree(uint64_t timestamp, addr_to_time_t* histogram)
// {
//   RDnode* new_tree = splay(0);
//   while (new_tree && new_tree->time < timestamp) {
//     RDnode* dead_node = new_tree;
//     new_tree = new_tree->right;
//     if (new_tree->left)
//       new_tree = new_tree->splay(0);
//     histogram->erase(dead_node->address);
//     delete dead_node;
//   }
//   return new_tree;
// }


// tree_dist() returns the number of nodes in a splay tree whose timestamp is
// larger than a given value.
uint64_t RDnode::tree_dist(uint64_t timestamp)
{
  RDnode* node = this;
  uint64_t num_larger = 0;
  while (true) {
    if (timestamp > node->time) {
      node = node->right;
    }
    else
      if (timestamp < node->time) {
        num_larger++;
        if (node->right != nullptr)
          num_larger += node->right->weight;
        node = node->left;
      }
      else {
        if (node->right != nullptr)
          num_larger += node->right->weight;
        return num_larger;
      }
  }
}


// Reinitialize an existing RDnode with a given address and timestamp.
void RDnode::initialize(uint64_t new_address, uint64_t new_time)
{
  address = new_address;
  time = new_time;
  weight = 1;
  left = nullptr;
  right = nullptr;
}


// Initialize a new RDnode with a dummy address and timestamp.
RDnode::RDnode()
{
  initialize(0, 0);
}


// Initialize a new RDnode with a given address and timestamp.
RDnode::RDnode(uint64_t address, uint64_t time)
{
  initialize(address, time);
}


// Define infinite distance.
const uint64_t infinite_distance = ~(uint64_t)0;


// A ReuseDistance encapsulates all the state needed for a
// reuse-distance calculation.
typedef std::unordered_map<uint64_t, uint64_t> addr_to_time_t;
class ReuseDistance {
private:
	addr_to_time_t last_access;   // Last access time of a given address
  uint64_t clock;           // Current time
  RDnode* dist_tree;        // Tree of reuse distances

public:
  // Initialize our various fields.
  ReuseDistance() {
    clock = 0;
    dist_tree = nullptr;
  }

  // Incorporate a new address into the reuse-distance histogram.
  inline uint64_t operator()(uint64_t address);

};


// Incorporate a new address into the reuse-distance histogram.
uint64_t ReuseDistance::operator()(uint64_t address)
{
  // Update the histogram.
  uint64_t distance = infinite_distance;
  addr_to_time_t::iterator prev_time_iter = last_access.find(address);
  RDnode* new_node = nullptr;
  if (prev_time_iter != last_access.end()) {
    // We've previously seen this address.
    uint64_t prev_time = prev_time_iter->second;
    distance = dist_tree->tree_dist(prev_time);
    dist_tree = dist_tree->remove(prev_time, &new_node);
  }


  // Update the tree and the map.
  if (new_node == nullptr)
    new_node = new RDnode(address, clock);
  else
    new_node->initialize(address, clock);
  if (__builtin_expect(dist_tree == nullptr, 0))
    // First insertion into the tree.
    dist_tree = new_node;
  else
    // All other tree insertions.
    dist_tree = dist_tree->insert(new_node);
  last_access[address] = clock;
  clock++;

  // If the tree and the map have grown too large, prune old addresses from
  // them.
//   if (last_access.size() > bf_max_reuse_distance)
//     dist_tree = dist_tree->prune_tree(clock - bf_max_reuse_distance, &last_access);
  return distance;
}  
