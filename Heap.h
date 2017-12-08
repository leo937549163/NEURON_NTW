/*
 * Heap.h
 *
 *  Created on: 2014-08-19
 *      Author: Mohammad Nazrul Ishlam Patoary
 */

#include <iostream>
#include <list>
#include <vector>
#include "SV.h"
using namespace std;

class Heap {
 public:
  long CurrentNum;
  SV* elist[1000000];  // Max limit of number of voxel !
  long elist_index[1000000];
  Heap();
  Heap(long MaxSize);
  ~Heap(void);
  SV* top(void);
  bool push(SV* Item);  // PUSH the Item to Heap
  SV* pop(void);        // POP and return Item from Heap
  long GetSize(void);   // Returns the number of nodes in the Heap
  void Display(void);
  vector<double> getTotal(int numberSpecies);
  long MAX_SIZE;                // The maximum number of elements
  void MoveUp(long Node);       // Shift Node up into place in PUSH
  void MoveDown(long Node);     // Shift Node down into place POP
  long ParentOf(long Node);     // Returns Parent location to MoveUp
  long LeftChildOf(long Node);  // Returns Left Child location to MoveDown
  void clear();
  void update(int targetID);
  int findTargetId(int targetID);
  bool LT(long double a, long double b);
};

Heap::Heap() {}

Heap::Heap(long MaxSize) : MAX_SIZE(MaxSize) { CurrentNum = 0; }

Heap::~Heap(void) {}

bool Heap::LT(long double a, long double b) {
  if (a < 0 || b < 0) {
    cout << "a: " << a << endl;
    cout << "b: " << b << endl;
    cout << "[ERROR][Heap.h] less than 0" << endl;
  }
  double epsilon = .0000000001;
  return (b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool Heap::push(SV* Item) {
  if (CurrentNum >= MAX_SIZE) {
    cout << "Limit exceed !" << endl;
    return false;
  }
  elist[CurrentNum] = Item;  // putting at the last position of Data[]
  elist_index[Item->ID] = CurrentNum;

  // cout << "Push"<< endl;
  // cout << "Item->ID:" << Item->ID << endl;
  // cout << "CurrentNum:" << CurrentNum << endl;
  // exit(-1);
  MoveUp(CurrentNum++);
  return true;
}

void Heap::MoveUp(long Node) {
  long Current = Node;
  long Parent = ParentOf(Current);
  // getting the index of current and Parent node
  SV* Item = elist[Current];
  while (Current > 0) {
    // if (elist[Parent]->pqTau > Item->pqTau){
    if (Item->pqTau < elist[Parent]->pqTau) {
      elist_index[elist[Parent]->ID] = Current;
      elist[Current] = elist[Parent];
      Current = Parent;
      Parent = ParentOf(Current);
    } else {
      break;
    }
  }  // end of while
  elist_index[Item->ID] = Current;
  elist[Current] = Item;
}

inline long Heap::ParentOf(long Node) {
  return (Node - 1) / 2;
}

SV* Heap /*<templ>*/ ::top(void) {
  SV* Temp = elist[0];
  return Temp;
}

SV* Heap /*<templ>*/ ::pop(void) {
  SV* Temp = elist[0];
  elist_index[CurrentNum] = 0;
  elist[0] = elist[--CurrentNum];  // Replace with the last Element
  MoveDown(0);
  return Temp;
}

void Heap::MoveDown(long Node) {
  long Current = Node;
  long Child = LeftChildOf(Current);
  // getting index of current and child
  SV* Item = elist[Current];          // Used to compare values
  while (Child < CurrentNum) {
    if (Child < (CurrentNum - 1)){
      if (elist[Child + 1]->pqTau < elist[Child]->pqTau) {
        Child++;
      }
    }
    if (elist[Child]->pqTau < Item->pqTau) {
      elist_index[elist[Child]->ID] = Current;
      elist[Current] = elist[Child];
      // Switch the Current node and the Child node
      Current = Child;
      Child = LeftChildOf(Current);
    } else {
      break;
    }
  }
  elist_index[Item->ID] = Current;
  elist[Current] = Item;
}

inline long Heap::LeftChildOf(long Node) {
  // return (Node * 2) + 1;
  return (Node << 1) + 1;
}

inline long Heap::GetSize(void) { return CurrentNum; }

void Heap::Display(void) {
  cout << "The Current HEAP " << endl;
  if (CurrentNum > 0) {
    cout << "{  Root->";
    for (int i = 0; i < CurrentNum; i++) {
      cout << " e:(";
      cout << elist[i]->getTau();
      cout << "," << elist[i]->pqTau;
      cout << "," << elist[i]->s;
      cout << "," << elist[i]->r;
      cout << "," << elist[i]->Lt;
      cout << "," << elist[i]->Gt;
      cout << ")" << "  ";
    }
    cout << " } " << endl;
  } else {
    cout << " Heap is EMPTY !" << endl;
  }
}  //

void Heap::update(int targetID) {
  int index = findTargetId(targetID);
  if (index < 0) {
    exit(-1);
  }
  MoveUp(index);
  MoveDown(index);
}

int Heap::findTargetId(int targetID) {
  return elist_index[targetID];
  // int index = -1;
  // for (int i = 0; i < CurrentNum; i++) {
  //   if (elist[i]->ID == targetID) {
  //     index = i;
  //
  //     // cout << (index==elist_index[targetID]);
  //
  //     // cout << "findTargetId"<< endl;
  //     // cout << "index:" << index << endl;
  //     // cout << "elist_index:" << elist_index[targetID] << endl;
  //     // exit(-1);
  //     return index;
  //   }
  // }
  // return index;
}

vector<double> Heap::getTotal(int numberSpecies) {
  vector<double> ss;
  if (CurrentNum > 0) {
    for (int s = 0; s < numberSpecies; s++) {
      double sum = 0;
      for (int i = 0; i < CurrentNum; i++) {
        sum = sum + elist[i]->Species[s]->numberOfMolecules;
      }
      ss.push_back(sum);
    }
  } else {
    cout << " Heap is EMPTY !" << endl;
  }
  return ss;
}

void Heap::clear() {
  for (int i = 0; i < CurrentNum; i++) {
    elist[i] = 0;
  }
  CurrentNum = 0;
}
