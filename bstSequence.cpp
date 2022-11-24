/*
 * Name: Lukas Bostick
 * Date Submitted: Nov 14, 2022
 * Lab Section: 003
 * Assignment Name: Lab 08
 */

#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <assert.h>

#include "bstSequence.h"

using namespace std;

void fix_size(Node *T)
{
  T->size = 1;
  if (T->left) T->size += T->left->size;
  if (T->right) T->size += T->right->size;
}

// insert value v at rank r
Node *insert(Node *T, int v, int r) {
  if (T == nullptr) return new Node(v);

  int tempRank = 0;
  if (T->left) { tempRank = T->left->size; }

  //if current rank is greater than or equal to target r, traverse down the left subtree.
  if (tempRank >= r) {
    T->left = insert(T->left, v, r);
  }
  //else traverse down the right tree, updating the rank. 
  else {
    T->right = insert(T->right, v, r - tempRank - 1);
  }

  fix_size(T);
  return T;
}

// returns a vector of key values corresponding to the inorder traversal of T (i.e., the contents of T in sorted order)
vector<int> inorder_traversal(Node *T)
{
  vector<int> inorder;
  vector<int> rhs;
  if (T == nullptr) return inorder;
  inorder=inorder_traversal(T->left);
  inorder.push_back(T->key);
  rhs=inorder_traversal(T->right);
  inorder.insert(inorder.end(), rhs.begin(), rhs.end());
  return inorder;
}

// return pointer to node of rank r (with r'th largest key; e.g. r=0 is the minimum)
Node *select(Node *T, int r)
{
  assert(T!=nullptr && r>=0 && r<T->size);

  int rank_of_root = T->left ? T->left->size : 0;
  if (r == rank_of_root) return T;
  if (r < rank_of_root) return select(T->left, r);
  else return select(T->right, r - rank_of_root - 1);
}

// Split tree T on rank r into tree L (containing ranks < r) and 
// a tree R (containing ranks >= r)
void split(Node *T, int r, Node **L, Node **R) {
  //base case
  if(T == nullptr) {
    *L = nullptr;
    *R = nullptr;
    return;
  }

  int tempRank = 0;
  if (T->left) { tempRank = T->left->size; }


  //If current rank is greater than rank r, that node and the right subtree will be included in tree R
  if (tempRank >= r) {
    //recursively call split down the path
    split(T->left, r, L, &T->left);
    //places the left child of current node T as the cell that *R points to, 
    //then moves pointer *R up a level to point to T 
    //(the top of the tree being rebuilt recursively).
    *R = T;
    //update size
    fix_size(*R);
  }
  //else that node and the left subtree will be included in tree L
  else {
    split(T->right, r - tempRank - 1, &T->right, R);
    *L = T;
    fix_size(*L);
  }
}

// insert value v at rank r
Node *insert_random(Node *T, int v, int r)
{
  // If (v,r) is the Nth node inserted into T, then:
  // with probability 1/N, insert (v,r) at the root of T
  // otherwise, insert_random (v,r) recursively left or right of the root of T
  if (T == nullptr) return new Node(v);
  
  double chance = rand() % (T->size);

  int tempRank = 0;
  if (T->left) { tempRank = T->left->size; }

  if (chance == 0) { 
    //create new empty trees L and R
    Node* L;
    Node* R;
    //set node containing k as new root.
    Node * newNode = new Node(v);

    split(T, r, &L, &R);
    //assign L and R to the appropriate subtrees of the new root (v, r)
    newNode->left = L;
    newNode->right = R;

    fix_size(newNode);
    return newNode;
  }
  else if (tempRank >= r) {
    T->left = insert_random(T->left, v, r);
    fix_size(T);
    return T;
  }
  else {
    T->right = insert_random(T->right, v, r - tempRank - 1);
    fix_size(T);
    return T;
  }
}

// Returns true if team x defeated team y
// Leave this function unmodified
bool did_x_beat_y(int x, int y)
{
  assert (x != y);
  if (x > y) return !did_x_beat_y(y,x);
  unsigned long long lx = x;
  unsigned long long ly = y;
  return ((17 + 8321813 * lx + 1861 * ly) % 1299827) % 2 == 0;
}

// Return a BST containing a valid ordering of n teams
Node *order_n_teams(int n)
{
  Node *T = nullptr;

  // start by inserting the first team
  T = insert_random(T, 0, 0);

  // now insert the other teams...
  for (int i=1; i<n; i++) {
    // insert team i so the sequence encoded by the BST remains valid
    if (did_x_beat_y(i, select(T,0)->key)) // can we insert at beginning?
      T = insert_random(T, i, 0);
    else if (did_x_beat_y(select(T,T->size-1)->key, i)) // can we insert at end?
	    T = insert_random(T, i, T->size);
    else {
      int low_rank = 0;
      int high_rank = T->size - 1;
      int mid_rank;
      
      //continues looping until the parameters are one position away from each other
      //(ensuring that the loop stops at a win, loss conjunction)
      while (high_rank - low_rank != 1) {
        mid_rank = (low_rank + high_rank) / 2;
        //if team to insert beat the team at mid_rank, adjust value of high_rank
        if (did_x_beat_y(i, select(T,mid_rank)->key)) {
          high_rank = mid_rank;
        }
        //otherwise, the team to insert lost to the team at mid_rank and low_rank must be adjusted 
        else {
          low_rank = mid_rank;
        }
      } 
      //now that the conjunction of loss/win is found, insert the team into the tree
      //before the team at high_rank
      T = insert_random(T, i, high_rank);
    }
  }
  return T;
}

void printVector(vector<int> v)
{
    for (int i=0; i<v.size(); i++)
    {
        cout << v[i] << endl;
    }
}

/*
int main(void)
{
  vector<int> inorder;
  Node *T = nullptr;

  // test insert at beginning
  for (int i=0; i<5; i++)
    T = insert(T, i+1, 0);
  cout << "Tree should contain 5 4 3 2 1:\n";
  inorder=inorder_traversal(T);
  printVector(inorder);

  // test insert at end
  for (int i=5; i<10; i++)
    T = insert(T, i+1, T->size);
  cout << "Tree should contain 5 4 3 2 1 6 7 8 9 10:\n";
  inorder=inorder_traversal(T);
  printVector(inorder);
  
  // test insert at middle
  for (int i=10; i<15; i++)
    T = insert(T, i+1, T->size/2);
  cout << "Tree should contain 5 4 3 2 1 12 14 15 13 11 6 7 8 9 10:\n";
  inorder=inorder_traversal(T);
  printVector(inorder);
    
  // once insert is working, the next step is to build the
  // insert_random function -- to test this, just change
  // calls to insert above to insert_random.

  int N = 100000; // this should run quickly even for very large N!
  Node *S = order_n_teams(N);
  if (S == nullptr || S->size != N)
    cout << "Size of tree returned by order_n_teams is wrong\n";
  else {
    cout << "Team ordering:\n";
    //    inorder=inorder_traversal(S);
    //    printVector(inorder);
    for (int i=0; i<N-1; i++) {
      Node *x = select(S, i);
      Node *y = select(S, i+1);
      if (!did_x_beat_y(x->key, y->key)) {
        cout << "Invalid sequence: team " << x->key << " (position " << i <<
	              ") lost to team " << y->key << " (position " << i+1 << ")\n";
      }
    }
  }  
  
  return 0;
}
*/

