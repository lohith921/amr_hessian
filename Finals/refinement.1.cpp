/* This program takes input the root of .node and .ele files  */

#include <map>
#include "mmap_lib_p.h"
#include "structures.h"

static int num_ele  =  0; // to keep track of total no of elements
static int num_nodes = 0; // to keep track of total no of nodes

/*************** Insert elements *****************************************************/
void insert_element(struct element *head_ele, struct element *t){
	struct element *temp;
	temp = new element();
	if(head_ele->next == NULL)
		head_ele->next = t;
	else{
		temp = head_ele->next;
		while(temp->next != NULL)
			temp = temp->next;
		temp->next = t;
	}
	num_ele++;
}

/***********************************************************************/
void print_element(struct element *e){
	printf("#: %d, Nodes %d, %d, %d\n",e->ele_no,e->nodes[0],e->nodes[1],e->nodes[2]);
}

/****************** Processes each element ****************************************************/
//void process(struct element *el, struct node_map *n_map, struct edge_map* emap, int tol){
void process(struct element *el, struct node_map *n_map, struct edge_map* emap){
     struct element *temp1, *temp2;
	 int i, j, eleno, vertices[3]; 
     temp1 = el;
     if (el != NULL){
		REAL *vertexA,*vertexB,*vertexC;
		REAL *midAB, *midBC, *midCA, length;
		vertexA = new REAL[2];
		vertexB = new REAL[2];
		vertexC = new REAL[2];
		struct element *nw, *nw2, *nw3;
		if (n_map != NULL){
			vertexA = map_getnode(n_map,temp1->nodes[0]);
			vertexB = map_getnode(n_map,temp1->nodes[1]);
			vertexC = map_getnode(n_map,temp1->nodes[2]);
		}
		else{
           std::cout << "The structures are not properly initiated" << std::endl;
           exit(0);
		}
		length = calc_length(vertexA, vertexB);
		//std::cout << "Length AB for " << temp1->ele_no << " is " << length << std::endl;
		if ( length > 1.414){
			midAB  =  new REAL[2]; 
			midAB  =  compute_mid(vertexA, vertexB);
			num_nodes++;
			map_setnode(n_map, num_nodes, midAB);
			nw = new element();
			nw->ele_no  =  ++num_ele;
			nw->nodes[0]  =  num_nodes;
			nw->nodes[1]  =  temp1->nodes[1];
			nw->nodes[2]  =  temp1->nodes[2];
			set_nedgemap(emap, nw->nodes);
			temp1->nodes[1]  =  num_nodes;
			set_nedgemap(emap, temp1->nodes);
			nw->next  =  temp1->next;
			temp1->next  =  nw; //new triangles A,mid,C and mid,B,C are formed
			//temp1 = nw;
			int ngb = edgemap_getnodes(emap, temp1->nodes[1], temp1->nodes[0]);
			//std::cout << "Neighboring node is " << ngb << std::endl;	        	
			if( ngb != -1){
				temp2 = find_element(el, temp1->nodes[0], temp1->nodes[1], ngb);
				if(temp2 != NULL){
					//std::cout << " Neighbor tri found is " << temp2->ele_no << std::endl;
					nw2 = new element();
					nw2->ele_no  =  ++num_ele;
					nw2->nodes[0]  = temp1->nodes[0];
					nw2->nodes[1]  = num_nodes;
					nw2->nodes[2]  = ngb;
					set_nedgemap(emap, nw2->nodes);
					temp2->nodes[0] = num_nodes;
					temp2->nodes[1] = temp1->nodes[1];
					temp2->nodes[2] = ngb;
					set_nedgemap(emap, temp2->nodes);
					nw2->next = temp2->next;
					temp2->next = nw2;
					//std::cout << "Inserted new elements " << nw->ele_no << " " << nw2->ele_no << std::endl;
					//temp2 = nw2;
				}	
				else{
					//std::cout << "No shared triangle found, returning " << std::endl;
					return;
				}
			}
			process(temp1, n_map, emap);
			return;
		}
			//std::cout << "Calling calc_length BC for element " << temp1->ele_no << std::endl;
		length = calc_length(vertexB, vertexC);
		//std::cout << "Length BC for " << temp1->ele_no << " is " << length << std::endl;
		if (length > 1.414){
        // printf("90\n");
			midBC  =  new REAL[2];
			midBC  =  compute_mid(vertexB, vertexC);
			num_nodes++;
			map_setnode(n_map, num_nodes, midBC);	
			nw  =  new element();
			nw->ele_no  =  ++num_ele;
			nw->nodes[0]  =  temp1->nodes[0];
			nw->nodes[1]  =  num_nodes;
			nw->nodes[2]  =  temp1->nodes[2];
			set_nedgemap(emap, nw->nodes);
			temp1->nodes[2]  =  num_nodes;
			set_nedgemap(emap, temp1->nodes);
			nw->next  =  temp1->next;
			temp1->next  =  nw;//new triangles A,B,mid and mid,C,A are formed
			//temp1  =  nw;
			int ngb = edgemap_getnodes(emap, temp1->nodes[2], temp1->nodes[1]);
			//std::cout << "Node retrieved is " << ngb << std::endl;
			if (ngb != -1){
				temp2 = find_element(el, temp1->nodes[1], temp1->nodes[2], ngb);
				if (temp2 != NULL){
					nw2 = new element();
					nw2->ele_no  =  ++num_ele;
					nw2->nodes[0]  = temp1->nodes[1];
					nw2->nodes[1]  = ngb;
					nw2->nodes[2]  = num_nodes;
					set_nedgemap(emap,nw2->nodes);
					//nw3 = new element();
					//nw3->ele_no  =  ++num_ele;
					temp2->nodes[0]  = num_nodes;
					temp2->nodes[1]  = ngb;
					temp2->nodes[2]  = temp1->nodes[2];
					set_nedgemap(emap,temp2->nodes);
					nw2->next = temp2->next;
					temp2->next = nw2;
					//std::cout << "Inserted new elements " << nw->ele_no << " " << nw2->ele_no << std::endl;
					//temp2 = nw2;
				}
				else{
				//	std::cout << "No shared triangle found, returning " << std::endl;
				//exit(0);
					return;
				}
			}
			//processBC(temp1, n_map, emap, vertexB, vertexC);  		
			process(temp1,n_map, emap);
			return;
		}
    		//std::cout << "Calling calc_length CA for element " << temp1->ele_no << std::endl;
    	length  =  calc_length(vertexC, vertexA);
		//std::cout << "Length CA for " << temp1->ele_no << " is " << length << std::endl;
    	if(length > 1.414){ //   printf("91\n");
    		midCA  =  new REAL[2];
			midCA  =  compute_mid(vertexC,vertexA); 
			num_nodes++;
			map_setnode(n_map, num_nodes, midCA);
			nw  =  new element();
			nw->ele_no  =  ++num_ele;
			nw->nodes[0] = temp1->nodes[0];
			nw->nodes[1] = temp1->nodes[1];
			nw->nodes[2] = num_nodes;
			set_nedgemap(emap,nw->nodes);
			temp1->nodes[0] = num_nodes;
			set_nedgemap(emap, temp1->nodes);
			nw->next = temp1->next;
			temp1->next = nw;//new triangles A,B,mid and mid,B,C are formed
			//temp1 = nw;
			int ngb = edgemap_getnodes(emap, temp1->nodes[0], temp1->nodes[2]);
			if (ngb != -1){
				temp2 = find_element(el, temp1->nodes[2], temp1->nodes[0], ngb);
				if (temp2 != NULL){
					nw2 = new element();
					nw2->ele_no  =  ++num_ele;
					nw2->nodes[0]  = temp1->nodes[2];
					nw2->nodes[1]  = ngb;
					nw2->nodes[2]  = num_nodes;
					set_nedgemap(emap,nw2->nodes);
    	
					temp2->nodes[0]  = ngb;
					temp2->nodes[1]  = temp1->nodes[0];
					temp2->nodes[2]  = num_nodes;
					set_nedgemap(emap,nw3->nodes);
					
					nw2->next = temp2->next;
					temp2->next = nw2;
					//temp2->next = nw3;
					//std::cout << "Inserted new elements " << nw->ele_no << " " << nw2->ele_no << std::endl;
					//temp1 = nw;
					//temp2 = nw2;
				}
				else{
					//std::cout << "No shared triangle found, returning " << std::endl;
					//exit(0);
					return;
				}
			}
			//processCA(temp1, n_map, emap, vertexC, vertexA);
		process(temp1,n_map, emap);
		return;
		}
    	//std::cout << "All lengths under limit" << std::endl;
		//temp1 = temp1->next;
		//temp1 = temp1->next;
		//process(temp1, n_map, emap, tol);
		return ;
	}
}

/************************************************************************
*************************************************************************/
int main(int argc, char **argv){   
	struct amrgeo *m;
	int i,j,nele,nnode,node_no;
	std::string filename;
	std::ifstream fnode, fele;
	//std::vector<std::string> tok;
	std::string tok[4];
	std::string line;
	REAL ideal_length;
	REAL vertices[3];
	int nbs[3];
	m = new amrgeo();
	struct element *head_ele = new element();
	struct node *head_node = new node();
	head_ele->ele_no = -1;
	struct edge_map *emap;
    std::cout << "Please input the root of the file names without any extension" << std::endl;
 	std::cin >> filename;
 	std::cout << "Entered file name is " << filename << std::endl;
 	m->elefilename = (filename + ".ele");
 	m->nodefilename = (filename + ".node");
    fele.open(m->elefilename, std::ios::in);
	fnode.open(m->nodefilename,std::ios::in); 
	if(!fele || !fnode){
		std::cout << "Files not found, exiting" << std::endl;
		exit(0);
	}
    //std::cout << "File opening successful " << std::endl;
    getline(fele, line);
	tokenize(line, tok);
	nele = std::stoi(tok[0]);
	//std::cout << "No of elements is " << nele << std::endl;
	emap = new edge_map();
	for (int i = 1; i <= nele; i++){
	     struct element *temp = new element();
		 getline(fele, line);
	   //  std::cout << "line read is " << line << std::endl;
	     tokenize(line, tok);
	    // std::cout << "Tokenized values " << tok[0] << " " << tok[1] << " " << tok[2] << " " << tok[3] << std::endl;
	     temp->ele_no = std::stoi(tok[0]);
	     temp->nodes[0] = std::stoi(tok[1]);
	     temp->nodes[1] = std::stoi(tok[2]);
	     temp->nodes[2] = std::stoi(tok[3]);
	     insert_element(head_ele, temp);
	     set_nedgemap(emap, temp->nodes);
	     //std::cout << "Edgemap is set" << std::endl;
	}
	fele.close();
	//std::cout << "Reading elements completed" << std::endl;
	//display_elements(head_ele);
	/* following code reads node file and creates node map. */
	
	char dmy1[3];
	char dmy2[1];
	/* creating a map of nodes */
	struct	node_map *nodemap  =  create_nodemap();
	getline(fnode,line);
	tokenize(line,tok);
	nnode = std::stoi(tok[0]);
	m->attribindi = std::stoi(tok[2]);
	m->boundindi = std::stoi(tok[3]);
	num_nodes  =  nnode;
	j = 0;
	for(int i = 1; i <= nnode; i++){
		getline(fnode,line);
 		tokenize(line,tok);
 		node_no = std::stoi(tok[0]);
 		vertices[0] = std::stod(tok[1]);
 		vertices[1] = std::stod(tok[2]);
 		if (m->boundindi)
 			vertices[2] = 1;	
 		//std::cout << "Tokenized values are" << node_no << " " << vertices[0] << " " << vertices[1] << std::endl;
 		map_setnode(nodemap, node_no, vertices);
	}
	//std::cout << "Reading nodes is done" << std::endl;
	fnode.close(); 
	//std::cout << "Displaying edge map" << std::endl;
	//display_edgemap(emap);
	struct element *temp = head_ele->next;
	std::cout << "Process is running" << std::endl;
	compute_quality(temp, nodemap);
	while (temp != NULL){
		//process(temp, nodemap, emap, ideal_length);
		//process(temp, nodemap, emap);
		std::cout << "Quality is" << temp->qual << std::endl;
		temp = temp->next;
	}
	head_ele = sort_list(head_ele);
	write_elements(filename, head_ele, num_ele);
	write_nodes(filename, nodemap, num_nodes);
	std::cout << "Writing to files successful" << std::endl;
	return 0;
}  
