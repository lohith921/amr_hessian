/*
 map_lib
 A simple associative-array library for C
*/

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <iostream>
#include "mmap_lib_p.h"

#ifdef SINGLE
#define REAL float
#else
#define REAL double
#endif

// The following function creates a map for nodes.
struct node_map *create_nodemap() {
   struct node_map *m;
   m = (struct node_map *)malloc(sizeof(struct node_map));
   if(!m)
      return NULL;
   m->point = 0;
   m->xandy = (REAL *)malloc(sizeof(REAL));
   m->xandy[0] = 0.0;
   m->xandy[1] = 0.0;
   m->nxt = NULL;
   return m;
}

// The following function frees node map.
void map_freenode(struct node_map *map) {
   if(!map)
      return;
   struct node_map *m, *mp;
   m = map;
   while(m != NULL) {
      mp = m;
      m = m->nxt;
      free(mp);
   }
}

// The following function is used to set node entries in node map.
void map_setnode(struct node_map *m,int pt,REAL coords[]) {
   struct node_map *map;
   // First time inserting
    if(m->point == 0){
    m->point = pt;
    m->xandy[0] = coords[0];
    m->xandy[1] = coords[1];
    m->nxt = NULL;
    return;
    }
    // #1 if pt matches any of the points already present
    for(map = m;;map = map->nxt) {
      if(pt == map->point){
          map->xandy[0] = coords[0]; 
	  map->xandy[1] = coords[1];           
          return;
      }// #1 ends here
      // #2 inserting at the last element
      if(map->nxt == NULL) {
          map->nxt = create_nodemap();
          //map->nxt=map_create();
          if(!map->nxt)
          	return;
	  map = map->nxt;
	  map->point = pt;
	  map->xandy[0] = coords[0];
	  map->xandy[1] = coords[1];
	  map->nxt = NULL;
          return;
      }// #2 ends here
     }
}

// The following function prints the node map.
void display_map(struct node_map *first){
 std::cout << "Displaying the Node map"  << std::endl;
 struct node_map *temp;
 temp = first;
 if(temp != NULL){
  while(temp != NULL){
   std::cout << "Node : " << temp->point << " x= " << temp->xandy[0] << " y= " << temp->xandy[1] << std::endl;
   temp = temp->nxt;
  }
 }
 else{
  std::cout << "Invalid map\n" << std::endl;
  return; 
 }
}

// The following function is used to get a node from the node map.
REAL* map_getnode(struct node_map *m,int pt) {
   struct node_map *map;
   for(map = m;map != NULL;map = map->nxt) {
     if(map->point == pt) //!strcasecmp(name,map->name)) {
       return map->xandy;
   }
   return NULL;
}

// The following function free's edge map.
void free_edgemap(struct edge_map *map) {
   if(!map)
      return;
   struct edge_map *m, *mp;
   m = map;
   while(m!=NULL) {
      mp = m;
      m = m->nxt;
      free(mp);
   }
}

// The following function is used to set node entries in edge map.
void set_edgemap(struct edge_map *m,int k1, int k2,int node) {
	//std::cout << "Keys are " << k1 << " " << k2 << " value is " << node << std::endl;
   // First time inserting
    if (m->node_num == -1){ //&& (m->node_nums[1] == -1)){
    	m->key_pair.first = k1;
    	m->key_pair.second = k2;
    	m->node_num = node;
    	m->nxt = NULL;
    	return;
    }
    struct edge_map *map;
    // #1 if key1 & key2 matches any of the edges already present
    for(map = m;map->nxt != NULL;map = map->nxt) {
      if  ((k1 == map->key_pair.first) && (k2 == map->key_pair.second)){ //|| ((k1 == map->key2) && (k2 == map->key1)) ){
      		map->node_num = node;       
          	return;
        }
      }
      // #1 ends here
      // #2 inserting at the last element
      if (map->nxt == NULL) {
          map->nxt = new edge_map();
         if(!map->nxt)
            return;
	 map = map->nxt;
	 map->key_pair.first = k1;  map->key_pair.second = k2;
	 map->node_num = node;
	 map->nxt = NULL;
         return;
      }// #2 ends here
}

void set_nedgemap(struct edge_map *m, int nodes_list[3]){
	//std::cout << "Set edge map is called with " << nodes_list[0] << " " << nodes_list[1] << " " << nodes_list[2] << std::endl;
	for(int i = 0;i < 3;i++)
		set_edgemap(m,nodes_list[i],nodes_list[(i+1)%3],nodes_list[(i+2)%3]);
}

void display_edgemap(struct edge_map *first){
	std::cout << "Displaying the edge map\n" << std::endl;
	struct edge_map *temp;
	temp = first;
	if (temp != NULL){
		while(temp != NULL){
		std::cout << temp->key_pair.first << "," << temp->key_pair.second << " Node is: " << temp->node_num <<std::endl;
		temp = temp->nxt;
		}
	}
	else{
		printf("Invalid map\n");
		return; 
	}
}
// The following function is used to get nodes from the edge map.
//std::vector<int> edgemap_getnodes(struct edge_map *m,int k1, int k2) {
int edgemap_getnodes(struct edge_map *m, int k1, int k2){   
   struct edge_map *map;
   for(map = m;map != NULL;map = map->nxt) {
     if ((map->key_pair.first == k1) && (map->key_pair.second == k2)) 
       return map->node_num;
   }
   return -1;
}

/*******************************************************************
// The following function creates a map for neighbors.
struct map_neighbor *create_neighbormap() {
   struct map_neighbor *m;
   m = (struct map_neighbor *)malloc(sizeof(struct node_map));
   if(!m)
      return NULL;
   m->t_num = 0;
   m->nbs = (int *)malloc(sizeof(int));
   m->nbs[0] = 0;
   m->nbs[1] = 0;
   m->nbs[2] = 0;
   m->nxt = NULL;
   return m;
}

// The following function is used to get neibhors from the neighbor map.
int* map_getneighbors(struct map_neighbor *m,int t) {
   struct map_neighbor *map;
   for(map = m;map != NULL;map = map->nxt) {
     if(map->t_num == t){ //!strcasecmp(name,map->name)) {
       return map->nbs;
     }
   }
   return NULL;
}

// The following function is used to set neighbor entries in neighbor map.
void map_setneighbor(struct map_neighbor *m,int t,int neighbs[]) {
   struct map_neighbor *map;
   // First time inserting
    if(m->t_num == 0){
    m->t_num = t;
    m->nbs[0] = neighbs[0];
    m->nbs[1] = neighbs[1];
    m->nbs[2] = neighbs[2];
    m->nxt = NULL;
    return;
    }
    // #1 if pt matches any of the points already present
    for(map = m;;map = map->nxt) {
      if(t == map->t_num){
          map->nbs[0] = neighbs[0]; 
	  map->nbs[1] = neighbs[1];
	  map->nbs[2] = neighbs[2];           
          return;
        }// #1 ends here
      // #2 inserting at the last element
      if(map->nxt == NULL) {
          map->nxt = create_neighbormap();
         //map->nxt=map_create();
         if(!map->nxt){
            return;
         }
	 map = map->nxt;
	 map->t_num = t;
	 map->nbs[0] = neighbs[0];
	 map->nbs[1] = neighbs[1];
	 map->nbs[2] = neighbs[2];
	 map->nxt = NULL;
         return;
      }// #2 ends here
   }
}

// The following function frees neighbor map.
void map_freeneighbor(struct map_neighbor *map) {
   if(!map)
      return;
   struct map_neighbor *m, *mp;
   m = map;
   while(m!=NULL) {
      mp = m;
      m = m->nxt;
      free(mp);
   }
}
/*em = (struct edge_map *)malloc(sizeof(struct edge_map));
   	if(!em)
      		return NULL;
	em->key1 = 0;
	em->key2 = 0;
   	em->nxt = NULL;
   	 //std::vector<int> dummy;
   //dummy = {0,0};
   /* The following function creates an edgemap.
struct edge_map* create_edgemap(){
	struct edge_map *em = new edge_map();
   	return em;
}
//m->key1 = k1;
    	//m->key2 = k2;
    	//m->node_nums.erase(m->node_nums.begin(), m->node_nums.end());
    	//map->node_nums.push_back(node); 
    //map->node_nums.erase(map->node_nums.begin(), map->node_nums.end());
	 //map->node_nums.push_back(node);
	 //|| ((map->key1 == k2) && (map->key2 == k1)) )
**************************************************************/
