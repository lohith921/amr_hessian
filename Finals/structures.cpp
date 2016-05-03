/* This is used to store all structure definitions */
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include "structures.h"
#include "mmap_lib_p.h"

/*************** Quality function ***************************/
void compute_quality(struct element* el, struct node_map* nmap){
	std::cout << "compute quality is called" << std::endl;	
	struct element* temp;
	REAL l1, l2, l3, A, p = 0.0, q;
	REAL *vertexA, *vertexB, *vertexC;
	vertexA = new REAL[2];
	vertexB = new REAL[2];
	vertexC = new REAL[2];
	temp = el;
	if(nmap != NULL){
		while (temp != NULL){
			vertexA = map_getnode(nmap, temp->nodes[0]);
			vertexB = map_getnode(nmap, temp->nodes[1]);
			vertexC = map_getnode(nmap, temp->nodes[2]);
			l1 = pow((vertexA[0] - vertexB[0]),2) + pow((vertexA[1] - vertexB[1]),2);
			l1 = sqrt(l1);
			l2 = pow((vertexB[0] - vertexC[0]),2) + pow((vertexB[1] - vertexC[1]),2);
			l2 = sqrt(l2);
			l3 = pow((vertexC[0] - vertexA[0]),2) + pow((vertexC[1] - vertexA[1]),2);
			l3 = sqrt(l3);
			p = (l1 + l2 + l3)/2;
			A = sqrt(p*(p-l1)*(p-l2)*(p-l3));
			q = 6.928203*A/(pow(l1,2) + pow(l2,2) + pow(l3,2));
			temp->qual = q;
			temp = temp->next;
		}
	}
	else{
         std::cout << "The structures are not properly initiated" << std::endl;
         exit(0);
	}
	
}
/**************String tokenizer for reading files*********************/
void tokenize(std::string str, std::string *token_v){
	//std::cout << "Tokenizer called " << std::endl;
	size_t start = str.find_first_not_of(" "), end = start;
	int i = 0;
	while (start != std::string::npos){
        // Find next occurence of delimiter
        	end = str.find(" ", start);
        // Push back the token found into vector
                //token_v.push_back(str.substr(start, end-start));
                token_v[i] = str.substr(start, end-start);
                i++;
        // Skip all occurences of the delimiter to find new start
        	start = str.find_first_not_of(" ", end);
    	}
}
/*******************************************************************/
REAL * compute_mid(REAL *vA, REAL *vB){
	REAL *nw;
	nw = new REAL[2];
	nw[0] = (vA[0]+vB[0]) / 2;
	nw[1] = (vA[1]+vB[1]) / 2;
	return nw;
}

/*******************************************************************/
REAL calc_length(REAL *A, REAL *B){
	REAL d1,d2,d;
	d1 = A[0]-B[0];
	d2 = A[1]-B[1];
	REAL temp = (d1*d1) + (d2*d2);
	d = std::sqrt(temp);
	return d;
}
/**********************************************************************/
void write_elements(std::string rt, struct element *el, int num_ele){
	std::ofstream fp; 
	fp.open((rt + ".r.ele"), std::ios::out);
	struct element *temp;
	temp = el->next;
	fp << num_ele << " " << 3 << std::endl; 
	while(temp !=  NULL){
		fp << temp->ele_no << " " << temp->nodes[0] << " " << temp->nodes[1] << " " << temp->nodes[2] << std::endl;
     		temp = temp->next;
	}
	fp.close(); 
}
/***********************************************************************/
void write_nodes(std::string rt,struct node_map *m, int num_nodes){
	std::ofstream fp;
	struct node_map *temp = m;
	fp.open((rt + ".r.node"),std::ios::out); 
	int i;
	fp << num_nodes << " " << 2 << std::endl; 
	while(temp!=  NULL){
	 	fp << temp->point << " " << temp->xandy[0] << " " << temp->xandy[1] << std::endl;
		temp = temp->nxt;  
	}
	fp.close();
}

/*********************************************************************/
void copy_nodes(int src[3],int dest[3]){
 int i;
 for(i = 0;i<3;i++){
  dest[i] = src[i];
 }
}

/*********************************************************************/
// This is very basic sorting, it requires improvement later
struct element *sort_list(struct element *head) {
    struct element *tmpPtr  =  head->next;
    struct element *tmpNxt  =  tmpPtr->next;
    int tmp;
    int tmp_nodes[3];
    while(tmpNxt !=   NULL){
           while(tmpNxt !=   tmpPtr){
                    if(tmpNxt->ele_no < tmpPtr->ele_no){
                            tmp  =  tmpPtr->ele_no;
                            copy_nodes(tmpPtr->nodes,tmp_nodes);
                           // tmp_nodes = tmpPtr->nodes;
                            tmpPtr->ele_no  =  tmpNxt->ele_no;
			    copy_nodes(tmpNxt->nodes,tmpPtr->nodes);
			    //tmpPtr->nodes = tmpNxt->nodes;
                            tmpNxt->ele_no  =  tmp;
			    copy_nodes(tmp_nodes,tmpNxt->nodes);
			   // tmpNxt->nodes = tmp_nodes;
                    }
                    tmpPtr  =  tmpPtr->next;
            }
            tmpPtr  =  head;
            tmpNxt  =  tmpNxt->next;
    }
         return tmpPtr ; // Place holder
} 

/*************************************************************************/
struct element* find_element(struct element *el, int n1, int n2, int n3){
	//std::cout << "Find element is called for nodes " << n1 << " " << n2 << " " << n3 << std::endl;
	struct element *temp;
	int nd1, nd2, nd3, c = 0;
	if(el->ele_no == -1)// for head element
                 temp = el->next;
     	else
                 temp = el;
     	//pv = temp;
	while(temp != NULL){
		c = 0;
		nd1 = temp->nodes[0]; nd2 = temp->nodes[1]; nd3 = temp->nodes[2];
		//std::cout << " Checking against following elements " << nd1 << " " << nd2 << " " << nd3 << std::endl;
		if ((nd1 == n1) || (nd1 == n2) || (nd1 == n3))
			c++;
		if ((nd2 == n1) || (nd2 == n2) || (nd2 == n3))
			c++;
		if ((nd3 == n1) || (nd3 == n2) || (nd3 == n3))
			c++;
		//std::cout << "Computed c value is " << c << std::endl;	
		if (c == 3)
			return temp;
		else
			temp = temp->next;
	}
	if (temp == NULL)
		return NULL;
}

/************************** Removed code *************************/
/*if ((nd1 == n1) && (nd2 == n2) && (nd3 == n3))
			return temp;
		else if ((nd1 == n2) && (nd2 == n3) && (nd3 == n1))
			return temp;
		else if ((nd1 == n3) && (nd2 == n1) && (nd3 == n2))
			return temp;
		else
		 	//pv = temp;
			temp = temp->next;
			*/
			/*for(int i = 0;i < 3;i++){
				vertices[i] = new REAL[2];
				vertices[i] = map_getnode(nmap,temp->nodes[i]);
			}				
			for (int j = 0;j < 3;j++){
				for( int k = 0;k < 2;k++){
					L[j] += pow((vertices[j][k] - vertices[(j+1)%3][k]),2);
				}
				L[j] = sqrt(L[j]);
				p += L[j];
			}*/
/****************** Set edge map************************
void set_nedgemap(std::map< std::pair<int,int>, int> m, int nodes_list[3]){
	//std::cout << "Set edge map is called with " << nodes_list[0] << " " << nodes_list[1] << " " << nodes_list[2] << std::endl;
	for(int i = 0;i < 3;i++)
		set_edgemap(m,nodes_list[i],nodes_list[(i+1)%3],nodes_list[(i+2)%3]);
}*/
