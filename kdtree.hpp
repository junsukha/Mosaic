/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>
#include <cmath> //for pow
using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */
  //std::cout << Dim << std::endl;
  if(first[curDim] == second[curDim])
    return first < second;
  return first[curDim] < second[curDim];
    //return false;
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     */
     int distance1 = distance(target, currentBest);
     int distance2 = distance(target, potential);
     if (distance2 == distance1)
      return potential < currentBest;
     return distance2 < distance1;

}
/////////////////////////////////////////////////
template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    /**
     * @todo Implement this function!
     */
     points = newPoints;
     if (newPoints.size() <= 1)
      return;
     KDTreeHelp(0, points.size()-1, 0);
     root = treeConstruct(root, points, 0, points.size()-1);
}

template <int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::treeConstruct(KDTreeNode* subroot, vector<Point<Dim>> points, int start, int end){
  int med =  (start + end) / 2;
  subroot = new KDTreeNode(points[med]);

  if(start >= end){
    subroot->left = nullptr;
    subroot->right = nullptr;
    //return;
    return subroot;
  }

  else if(end-start == 1){
    subroot->left = nullptr;
    subroot->right = treeConstruct(subroot->right, points, med+1, end);
  }

  else{
    subroot->left = treeConstruct(subroot->left, points, start, med-1);
    subroot->right = treeConstruct(subroot->right, points, med+1, end);
  }

  return subroot;
}


template <int Dim>
void KDTree<Dim>::KDTreeHelp(int start, int end, int dim){
  if (start>=end)
    return;
  int med = (start + end) / 2;
  quickSelect(start, end, med, dim);
  KDTreeHelp(start, med-1, (dim+1) % Dim);
  KDTreeHelp(med+1, end, (dim+1) % Dim);
}

template <int Dim>
void KDTree<Dim>::quickSelect(int start, int end, int med, int dim)
{
  if (start >= end) return;

  int index = partition(start, end, (start + end) / 2, dim);
  if (index > med)
    return quickSelect(start, index-1, med, dim);
  if (index < med)
    return quickSelect(index+1, end, med, dim);
}

template <int Dim>
int KDTree<Dim>::partition(int start, int end, int pIndex, int dim){
  Point<Dim> pivot = points[pIndex];
  std::swap(points[pIndex], points[end]);
  int index = start;
  for (int i = start; i < end; i++) {
    if(smallerDimVal(points[i], pivot, dim)){
      std::swap(points[index], points[i]);
      index++;
    }
  }
  std::swap(points[end], points[index]);
  return index;
}

/* partition puedocode
template <int Dim>
int KDTree<Dim>::partition(list, left, right, pivotIndex){
     pivotValue := list[pivotIndex]
     swap list[pivotIndex] and list[right]  // Move pivot to end
     storeIndex := left
     for i from left to right-1
         if list[i] < pivotValue
             swap list[storeIndex] and list[i]
             increment storeIndex
     swap list[right] and list[storeIndex]  // Move pivot to its final place
     return storeIndex
}
*/
///////////////////////////////////////////



template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
   copy(this->root, other.root);
   size = other.size;
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   */
   if (this != NULL) clear(root);
   copy(this->root, rhs.root);
   size = rhs.size;

   return *this;
}

template <int Dim>
void KDTree<Dim>::copy(KDTreeNode * subroot, KDTreeNode * othersubroot) {
	subroot = new KDTreeNode();
	subroot->point = othersubroot->point;
	copy(subroot->left, othersubroot->left);
	copy(subroot->right, othersubroot->right);
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
   //clear(root);
   //clear(root);
   //cleantree(root);
}

template <int Dim>
void KDTree<Dim>::clear(KDTreeNode * node) {
  if (node)
  {
      clear(node->left);
      clear(node->right);
      delete node;
  }

  // if(subroot->left == nullptr && subroot->right == nullptr)
  //   free(subroot);
  // clear(subroot->left);
  // clear(subroot->right);


	// if (subroot == NULL) return;
	// if (subroot->left != NULL) clear(subroot->left);
	// if (subroot->right != NULL) clear(subroot->right);
	// delete subroot;
	// subroot = NULL;
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
    //query is the target point. Find a point that is the closest to the target point from
    //my "points" vector
    //Point<Dim> result;
    int result;
    //compare with the root


    // if(query == root)
    //   return root;


    // else if(query[0] < root[0])
    //   result = findNearestNeighbor(query, root->left);
    // else
    //   result = findNearestNeighbor(query, root->right);

    //start from dim = 0
    result = findNearestNeighbor(query, 0, points.size() - 1, 0);//points.size()
    return points[result];
    //return Point<Dim>();
}


template <int Dim>
int KDTree<Dim>::findNearestNeighbor(const Point<Dim> query, int left, int right, int dim) const{
  if(left >= right)
    //return right;
    return left;


  int currentBest;
  int med = (left + right)/2;
  // if(left >= right)
  //   return left;
  if(smallerDimVal(query, points[med], dim))//query < points[med]
    currentBest = findNearestNeighbor(query, left, med - 1, (dim + 1) % Dim);
  else
    currentBest = findNearestNeighbor(query, med + 1, right, (dim + 1) % Dim);


  //now start back traversal to check for other nodes within this radius
  //int radius =
  //shoudReplace(target, currentBest, potential)

  //back at parent node, compare the distance with curretBest and update it the closer one.
  if(shouldReplace(query, points[currentBest], points[med]))
    currentBest = med;

  //update the radius if changed
  int newDistance = distance(query, points[currentBest]);


  int splitDis = points[med][dim] - query[dim];
  //distance is using r^2 not r. So I need to squre the splitDis also.
  splitDis = pow(splitDis, 2);

  //check if splitting plane is within the radius
  int temp;
  if(splitDis <= newDistance){
    if(smallerDimVal(query, points[med], dim))
      //int temp = findNearestNeighbor(query, med + 1, right, (dim + 1) % Dim);
      temp = findNearestNeighbor(query, med + 1, right, (dim + 1) % Dim);
    else
      //int temp = findNearestNeighbor(query, left, med - 1, (dim + 1) % Dim);
      temp = findNearestNeighbor(query, left, med - 1, (dim + 1) % Dim);

    if(shouldReplace(query, points[currentBest], points[temp]))
      currentBest = temp;

    //return currentBest;
  }
  return currentBest;

}
// template <int Dim>
// Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim> query, KDTreeNode* subroot){
//
// }

//custom helper function
template <int Dim>
int KDTree<Dim>::distance(const Point<Dim> first, const Point<Dim> second) const{
  int result = 0;
  for(int i = 0; i < Dim; ++i){
    result += pow((first[i] - second[i]), 2);
  }
  return result;
}
