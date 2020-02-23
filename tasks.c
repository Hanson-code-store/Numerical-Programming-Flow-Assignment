/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : 913022
 *   Name        : Thomas Hanson
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <assert.h>
#include "tasks.h"

#include <unistd.h>

#define RHO 1
#define U 2
#define V 3
#define X 4
#define Y 5
#define FLUX_U 6
#define FLUX_V 7
#define SCORE 8

#define MIN_X -15.0 
#define MAX_X 85.0
#define MIN_Y -20.0
#define MAX_Y 20.0

#define THRESH_1 5 
#define THRESH_2 10
#define THRESH_3 15
#define THRESH_4 20
#define THRESH_5 25

#define MIN_X_FOR_FLUX 20.0
#define ARRAY_SIZE 10000
#define ENTRY_SIZE 1000
#define MULTI 1000.0
#define BIN_SIZE 20
#define TWO_D  2
#define ON_LINE 1.0
#define X_COORD 0
#define Y_COORD 1

#define MIN 1
#define MAX 2
#define TRUE 1
#define FALSE 0
#define ERROR -1



/*Type defining structs for all functions */
typedef struct bst bst_t;
typedef struct node node_t;
typedef struct data data_t;

struct bst {
    int num_elements;
	node_t* root;
	
};

struct node {
	
    node_t* left;
    node_t* right;
    data_t* data;
};

struct data{
	double rho;
	double u;
	double v;
	double x;
	double y;
	double flux_u;
	double flux_v;
	double score;
	double w;
	data_t* nextdata;
	
};


data_t* getData(FILE *ifp);
bst_t* newBST(void);
node_t* makeMesh(int );
void addToArray2D(bst_t*, data_t*);
void balancedInsert(bst_t* , data_t** , int, int, int);
void nodeInsert(bst_t*, data_t*, int);
node_t* findMinMax(bst_t*, int);
void sortArrayA(data_t**, int, int);
double dataCompare(data_t*, data_t*, int);
int getGridNum(data_t*, int, int);
data_t* getAvg(data_t*);
data_t* putInList(int, data_t*, data_t*, char);
double arraySearchFluxU(data_t**, int, double,FILE *ofp);
int binarySearch(data_t**, int, int, double,FILE *ofp);
int searchLinkedList(data_t* data_list,double target_flux_u, FILE *ofp);
int searchBST(bst_t* , double, FILE *ofp);
void vortex(bst_t*, int, int, int*);
data_t* bstT(bst_t*, int, int);
char* fgetstr(char*, int, FILE *stream);
void findLength(bst_t*, int*);
void printNode(node_t*, FILE *ofp);
void printDataS(data_t*, FILE *ofp);

double getMicrotime(void);
void freeBST(bst_t*,int);
void freeTree(node_t*,int);
void freeData(data_t*);

void maxfluxdiff(const char* flow_file)
{
	FILE *ifp = NULL;
	FILE *ofp = NULL;
    data_t* new_data = NULL;
    bst_t* bst_flux_u = NULL;
    bst_t* bst_flux_v = NULL;
    char bin[BIN_SIZE];
    

//Opening data file
    ifp = fopen(flow_file,"r");
    assert(ifp); 
    
//Creating Output file
    ofp = fopen("task1.csv","w");
    assert(ifp); 
    
//Creating new BSTs
    bst_flux_u =(bst_t*)newBST();
    bst_flux_v =(bst_t*)newBST();
   	
//Reading data from file into 2 BST structures
    new_data = getData(ifp);
    fgets(bin,BIN_SIZE,ifp);
    while(new_data){
    	
    	if((new_data->x)>MIN_X_FOR_FLUX){ 
    		nodeInsert(bst_flux_u,new_data,FLUX_U);
    		nodeInsert(bst_flux_v,new_data,FLUX_V);
        }else{
        	free(new_data);
        }
        
        new_data = getData(ifp);
    }
    
//Finding minimum and maximum values and printing out values    
    fprintf(ofp,"rho,u,v,x,y\n");
    printNode(findMinMax(bst_flux_u,MAX),ofp);
    printNode(findMinMax(bst_flux_u,MIN),ofp);
    printNode(findMinMax(bst_flux_v,MAX),ofp);
    printNode(findMinMax(bst_flux_v,MIN),ofp);
   
//Closing files and freeing data
    fclose(ifp);
    fclose(ofp);
    freeBST(bst_flux_u,TRUE);
    freeBST(bst_flux_v,FALSE);
	return;
}

void coarsegrid(const char* flow_file, int resolution)
{
	FILE *ifp = NULL;
	FILE *ofp = NULL;
    data_t* new_data = NULL;
    data_t* current_data = NULL;
    node_t* current_node = NULL;
    node_t* head_node = NULL;
    bst_t* mesh = NULL;
    data_t* avg_data_list = NULL;
    
    char bin[BIN_SIZE];
    int grid_coord [TWO_D];
    int i = 0;
    int j = 0;
    
    
    
//Opening data file
    ifp = fopen(flow_file,"r");
    assert(ifp); 
    
//Creating Output file
    ofp = fopen("task2.csv","w");
    assert(ofp); 
    
//Creating new mesh
    mesh = (bst_t*)newBST();
    mesh->root = makeMesh(resolution);
    
//Reading data from file into grid structure
	fgets(bin,BIN_SIZE,ifp);
    new_data = getData(ifp);
    
    while(new_data){
    	i = 0;
    	
//Finds the x and y coords of the new data   	
    	grid_coord[X_COORD] = getGridNum(new_data,X,resolution);
    	grid_coord[Y_COORD] = getGridNum(new_data,Y,resolution);
 
//Checks if data lies copletely within a grid, then travels to that location    	
    	if((grid_coord[X_COORD] != -1) && (grid_coord[Y_COORD] != -1)){
    		current_node = mesh->root;
    		for(i = 0; i < grid_coord[X_COORD]; i++){
    			current_node = current_node->right;
    		}
    		for(i = 0; i < grid_coord[Y_COORD]; i++){
    			current_node = current_node->left;
    		}
    		
//Checks if it is the first data in that grid cell, if it isn't, the data
//is added to a linked list 
    		if(current_node->data == NULL){	
    			current_node->data = new_data;
    		}else{
    			current_data = current_node->data;
    			while(current_data->nextdata != NULL){
    				current_data = current_data->nextdata;
    			}
    			current_data->nextdata = new_data;
    		}
    	}else{
    		free(new_data);
    	}
    	new_data = getData(ifp);			
    }
    
//Calculates data average and puts data into sorted linked list
	head_node = mesh->root;
	
    for(i = 0; i < resolution; i++){
    	current_node = head_node;
    	for(j = 0; j < resolution; j++){
    			new_data = getAvg(current_node->data);
    			avg_data_list = putInList(SCORE,avg_data_list,new_data,'d');
    			current_node = current_node->left;
    	}
    	head_node = head_node->right;
    }
    
//Printing averaged data in order
    fprintf(ofp,"rho,u,v,x,y,S\n");
    current_data = avg_data_list;
    while(current_data){
    	printDataS(current_data,ofp);
    	current_data = current_data->nextdata;
    }
    
    freeBST(mesh,FALSE);
    freeData(avg_data_list);
    
    return;
}

void searching(const char* flow_file)
{
	int array_len = 0;
	int i = 0;
	double target_flux_u = 0;
	double start = 0;
	double stop = 0;
   
    char bin[BIN_SIZE];
    
	FILE *ifp = NULL;
	FILE *ofp = NULL;
    data_t* array[ARRAY_SIZE] = {NULL};
    data_t*  new_data = NULL;
    data_t* data_list = NULL;
    bst_t* ballenced_bst = NULL;
    
   
 
//Opening data file
    ifp = fopen(flow_file,"r");
    assert(ifp); 
    
//Creating Output file
    ofp = fopen("task3.csv","w");
    assert(ofp);
    
//Creating new BST
    ballenced_bst = newBST();
      
//Puts data that meets requrement (data->y = 0) into an array    
    fgets(bin,BIN_SIZE,ifp);
    new_data = getData(ifp);
    
     while(new_data){
     	if(new_data->y == 0){	
     		array[i] = new_data;
     		i++;
     	}else{
     		free(new_data);
     	}
     	new_data = getData(ifp);
     }
     array_len = i - 1;
     
     sortArrayA(array,array_len,FLUX_U);
     
//Sorts array into a linked list    
     for(i = 0; i < array_len; i++){
    	data_list = putInList(FLUX_U,data_list,array[i],'a');
    }

//Puts data into a ballenced BST
    balancedInsert(ballenced_bst,array,0,array_len,FLUX_U);
    
    
//Finding value with 40% max flux   
	start = getMicrotime();
	target_flux_u = arraySearchFluxU(array, array_len, 0.4, ofp); 
	stop = getMicrotime();
    fprintf(ofp,"\n");
	printf("TASK 3 Array Linear Search:  %.2lf microseconds\n", stop - start);

//Performing a binary search on the data
	start = getMicrotime();
    binarySearch(array, 0, array_len, target_flux_u, ofp);
    stop = getMicrotime();
    fprintf(ofp,"\n");
    printf("TASK 3 Array Binary Search:  %.2lf microseconds\n", stop - start);
   
//Performing a linear search on the linked list
	start = getMicrotime();
	searchLinkedList(data_list,target_flux_u, ofp);
	stop = getMicrotime();
	fprintf(ofp,"\n");
    printf("TASK 3 List Linear Search:  %.2lf microseconds\n", stop - start);
	
//Performing a  search on the BST
	start = getMicrotime();
	searchBST(ballenced_bst, target_flux_u, ofp);
	stop = getMicrotime();
	fprintf(ofp,"\n");
    printf("TASK 3 BST Search:  %.2lf microseconds\n", stop - start);
    
    
    
    
    freeData(data_list);
    freeBST(ballenced_bst,FALSE);
}

void vortcalc(const char* flow_file)
{
    char bin[BIN_SIZE];
    int xy[TWO_D];
    int y = 0;
    int x = 0;
    int thresh_five = 0;
    int thresh_ten = 0;
    int thresh_fift = 0;
    int thresh_twenty = 0;
    int thresh_twentyf = 0;
    double w = 0;
    
    FILE *ifp = NULL;
	FILE *ofp = NULL;
    data_t* new_data = NULL;
    bst_t* twodem_bst = NULL;
    
//Opening data file
    ifp = fopen(flow_file,"r");
    assert(ifp); 
    
//Creating Output file
    ofp = fopen("task4.csv","w");
    assert(ofp);
    
//Creating new BST
    twodem_bst = newBST();
    
//Adds data to a 2D dynamic array using the bst structure.
//Moving right is analogous to traversal along the y axis and vice versa
    fgets(bin,BIN_SIZE,ifp);
    new_data = getData(ifp);
    
     while(new_data){
     	addToArray2D(twodem_bst, new_data);
     	new_data = getData(ifp);
     }
     
//Finds size of the 2D array     
     findLength(twodem_bst, xy);
     
     
//For every point in the array, calculates the w value and finds its threshold     
     for(y = 0; y < xy[Y_COORD]; y++){
     	 for(x = 0; x< xy[X_COORD]; x++){
     	 	 vortex(twodem_bst,x,y,xy);
     	 	 
     	 	 w = fabs(bstT(twodem_bst,x,y)->w);
     	 	 if(w < THRESH_1){
     	 	 	 thresh_five++;
     	 	 }else if(w < THRESH_2){
     	 	 	 thresh_ten++;
     	 	 }else if(w < THRESH_3){
     	 	 	 thresh_fift++;
     	 	 }else if(w < THRESH_4){
     	 	 	 thresh_twenty++;
     	 	 }else if(w < THRESH_5){
     	 	 	 thresh_twentyf++;
     	 	 }
     	 }
     }
     
//Prints threshold data     
     fprintf(ofp,"threshold,points\n");
     fprintf(ofp,"5,%d\n10,%d\n15,%d\n20,%d\n25,%d\n",thresh_five,thresh_ten,
     	 thresh_fift,thresh_twenty,thresh_twentyf);
     
     freeBST(twodem_bst,TRUE);
     
     return;
}
     

//Funtions for tasks.c

/*Reads in data from file stream and puts it into a data_t struct.
 *Returns data_t pointer to new data struct.
 *Requires file pointer as an input argument. 
 */
data_t* getData(FILE *ifp){
	char buffer[ENTRY_SIZE];
	char entry[ENTRY_SIZE];
	int i = 0;
	int j = 0;
	int k = 1;

//Allocates memory for new daya	
    data_t* new = (data_t*)malloc(sizeof(data_t));
    assert(new);
    new->rho = 0;
    new->u = 0;
    new->v = 0;
    new->x = 0;
    new->y = 0;
    new->flux_u = 0;
    new->flux_v = 0;
    new->score = 0;
    new->w = 0;
    new->nextdata = NULL;
    
//Checks if at end of file    
    if (fgetstr(buffer,sizeof(char)*ENTRY_SIZE,ifp) == NULL){
        free(new);
        return NULL;
    }
    
//Reading in data
    while(buffer[i] != '\0'){
    	
    	if(buffer[i] == ','){
    		entry[j] = '\0';
    		if(k == RHO){
    			new->rho = strtod(entry,NULL);
    		}else if(k == U){
    			new->u = strtod(entry,NULL);
    			new->flux_u = (new->rho)*(new->u);
    		}else if(k == V){
    			new->v = strtod(entry,NULL);
    			new->flux_v = (new->rho)*(new->v);
    		}else if(k == X){
    			new->x = strtod(entry,NULL);	
    		}
    		i++;
    		j = 0;
    		k++;
    	}else{
    		entry[j] = buffer[i];
    		i++;
    		j++;
    	}
    }
  
//adds last entry    
    entry[j] = buffer[i];
    new->y = strtod(entry,NULL);
    
    return new;
}


/*Creates a new BST. Function returns pointer to BST.
 */
bst_t* newBST(void){
	bst_t* new_bst;
	
//Allocates memory for BST and initializes values
    new_bst = (bst_t*)malloc(sizeof(bst_t));
	assert(new_bst);
	new_bst->root = NULL;
	new_bst->num_elements = 0;
	
	return new_bst;
}

/*Creats a new mesh of size resolution x resolution. 
 *Requires int resolution argument and returns pointer to root node of mesh. 
 */
node_t* makeMesh(int resolution){
	int i = 0;
	int j = 0;
	node_t* current = NULL;
	node_t* head = NULL;

//Allocates memory for root node	
	node_t* root = (node_t*)malloc(sizeof(node_t));
	assert(root);
	root->left = NULL;
	root->right = NULL;
	root->data = NULL;
	
	head = root;

//Creates a resoltion x resolution grid of node_t 	
	for(i = 0; i < resolution; i++){
		current = head;
		for(j = 0; j < resolution; j++){
			current->left = (node_t*)malloc(sizeof(node_t));
			assert(current->left);
			current = current->left;
			current->left = NULL;
			current->right = NULL;
			current->data = NULL;
		}
		head->right = (node_t*)malloc(sizeof(node_t));
		assert(head->right);
		head = head->right;
		head->left = NULL;
		head->right = NULL;
		head->data = NULL;
	}
	return root;
}


/*Adds data to "2D" BST.
 *Requres pointer to data_t and pointer to bst_t as arguments.
 */
void addToArray2D(bst_t* bst, data_t* data){
	node_t* parent;
	
	node_t* new = (node_t*)malloc(sizeof(node_t));
	assert(new);
	new->left = NULL;
	new->right = NULL;
	new->data = data;
	bst->num_elements += 1;

//Tests for first node_t 	
	parent = bst->root;
	if(!parent){
		bst->root = new;
		return;
	}
//Looks along the 'y' axis until its row or a new row is found
	while(new->data->y != parent->data->y){
		if(parent->right == NULL){
			parent->right = new;
		}else{
			parent = parent->right;
		}
	}
//Looks along the 'x' axis until its position is found
	while(new->data->x != parent->data->x)
		if(parent->left == NULL){
			parent->left = new;
		}else{
			parent = parent->left;
		}
	return;
}


/*Inputs midpoints of array into BST such that the BST remains ballenced.
 *Requires pointer to bst_t and pointer to data_t array as arguments.
 *Function addapted from Chitrarth Lav, 2 August 2017.
 */ 
void balancedInsert(bst_t* bst, data_t* array[ARRAY_SIZE], 
int min, int max, int sort_key){
	int mid = 0;
	
	bst->num_elements += 1;
	if(min<=max){
		mid = min + (max - min) / 2;
		nodeInsert(bst,array[mid],sort_key);
		balancedInsert(bst,array,min,mid - 1,sort_key);
		balancedInsert(bst,array,mid+1,max,sort_key);
	}
	return;
}
			

/*Inserts new node into BST. 
 * Requires bst_t pointer, data_t pointer and integer sort and 
 * dupelacate keys as arguments.
 */ 
void nodeInsert(bst_t* bst, data_t* data, int sort_key){
   
    double comp;
    node_t* parent = NULL;
    
//Creates and initializes new node  
    node_t* new_node = (node_t*)malloc(sizeof(node_t));
    assert(new_node);
    new_node->left = NULL;
    new_node->right = NULL;
    new_node->data = data;
    
    bst->num_elements += 1;
    parent = bst->root;
    
//Check for first node in BST    
    if(!parent){
        bst->root = new_node;
        return;
    }
     
//Compares sort key of new data to parent data. 
//Then traverses tree until empty spot is found.    
    while(TRUE){
        comp = dataCompare(new_node->data,parent->data,sort_key); 
        
        if(comp<0){
            if(parent->left == NULL){
                parent->left = new_node;
                return;
            }else{
                parent = parent->left;
           }
        }else{
            if(parent->right == NULL){
                parent->right = new_node;
                return;
            }else{
            	
                parent = parent->right;
            }   
        }
    }
    return;
}


/*Sorts elements into a linked list. Data is linked to within data->nextdata.
 *Requires an int sort key, data_t root data, data_t new data and a char to 
 *assert assending or decending sort order.
 *Returns root data_t.
 */
data_t* putInList(int sort_key, data_t* parent, data_t* new_data,char updown){
	int comp;
    
    if(!parent){
        return new_data;
    }else{
//Compares new data to parent data based on the sort key    	
    	comp =  dataCompare(new_data,parent,sort_key);
    	
//Assending case. If new data < parent data, puts parent as the next element 
//after new data and returns new data to the previous parent.
//Otherwise it tries the next data entry recursively
    	if(updown == 'a'){
    		if(comp<0){
    			new_data->nextdata = parent; 
            	return new_data;
            }else{
            	parent->nextdata = 
    			putInList(sort_key,parent->nextdata,new_data,updown);
    			return parent;
        }
//Decending case
    	}else if(updown == 'd'){
    		if(comp<0){
    			parent->nextdata = 
    			putInList(sort_key,parent->nextdata,new_data,updown);
    			return parent;
    		}else{
    			new_data->nextdata = parent; 
    			return new_data;
    		}
    	}	
    }
    return NULL;
}


/*Sorts array in assending order.
 *requres data_t array[], int length of array and int sort key as inputs.
 */
void sortArrayA(data_t* array[ARRAY_SIZE],int length, int sort_key){
	int i = 0;
	int j = 0;
	data_t* tmp = NULL;

	for (i = 0; i < length; i++){
		for (j = 0; j < length; j++){
//If new data is larger than old data, switch their position
			if(dataCompare(array[j],array[i],sort_key)>0){
				tmp = array[i];         
				array[i] = array[j];           
				array[j] = tmp;             
			}  
		}
	}
	return;
}


/*Compares 2 data_t structs based on a sort key
 *Requres 2 data_t pointers and an int sort key as arguments.
 *Returns a double such that if new data is less than old data, return is < 0.
 */
double dataCompare(data_t* new, data_t* parent, int sort_key){
	double comp = 0;
	switch(sort_key){
		case RHO:
			comp =  new->rho - parent->rho;
			break;
		case U:
			comp =  new->u - parent->u;
			break;
		case V:
			comp =  new->v - parent->v;
			break;
		case X:
			comp =  new->x - parent->x;
			break;
		case Y:
			comp =  new->y - parent->y;
			break;
		case FLUX_U:
			
//Multiplier used in this case as flux values are typically very small
			comp = MULTI*(new->flux_u) - MULTI*(parent->flux_u);
			break;
		case FLUX_V:
			comp =  MULTI*(new->flux_v) - MULTI*(parent->flux_v);
			break;
		case SCORE:
			comp =  new->score - parent->score;
			break;
	}
	return comp;
}
	

/*Traverses tree left or right to find min or max value.
 *Requires pointer to bst_t and int min or max key as input arguments.
 *Returns pointer to node_t that has min or max value
 */
node_t* findMinMax(bst_t* bst, int min_max){
	node_t* parent = bst->root;
	
//Traverses tree left or right until NULL pointer is found	
	if(min_max == MIN){
		while(parent->left != NULL){
			parent = parent->left;
		}
		return parent;
	}else if(min_max == MAX){
		while(parent->right != NULL){
			parent = parent->right;
		}
		return parent;
	
	}else{
		return NULL;
	}
}


/*Finds the grid coordinates of a data_t struct within a mesh grid
 *Requires data_t pointer, int that states what coordinate is required and 
 *int resolution of grid.
 *Returns int which is the x or y index of the data_t within the mesh grid
 */
int getGridNum(data_t* data,int coord,int resolution){
	double norm = 0;
	
//Xvalue is requested. Data is normalised between 0 and 1, then multiplied by
//the resolution
	if(coord == X){
		norm = (data->x - MIN_X)/(MAX_X - MIN_X);
		norm = norm * resolution;
		
//Yvalue is requested		
	}else if(coord == Y){
		norm = (double)(data->y - MIN_Y)/(double)(MAX_Y - MIN_Y);
		norm = norm * resolution;
		
	}
	
//Checks if data sits on grid line, if so return int ERROR	
	if(fmod(norm,ON_LINE) == 0){
		return ERROR;
	}else{
		return (int)floor(norm);
	}
	return ERROR;
}


/*Gets the average values of all data_t structs within a linked list and 
 *calulates a score for the average.
 *Requires data_t pointer to root of linked list.
 *Returns data_t pointer to data_t with average values.
 */
data_t* getAvg(data_t* current_data){
	double rho_t = 0;
	double u_t = 0;
	double v_t = 0;
	double x_t = 0;
	double y_t = 0;
	int length = 0;
	data_t* previous = NULL;
	

//Adds data_t values to total values, moves to next in linked list and
//records length of list
	while(current_data->nextdata != NULL){
		rho_t = rho_t + current_data->rho;
		u_t = u_t + current_data->u;
		v_t = v_t + current_data->v;
		x_t = x_t + current_data->x;
		y_t = y_t + current_data->y;
		length++;
		previous = current_data;
		current_data = current_data->nextdata;
		free(previous);
	}

//For last data_t struct, calculates averages and score, then returns its 
//data_t pointer
	length = (double)(length + 1);
	current_data->rho = (rho_t + current_data->rho)/(length);
	current_data->u = (u_t + current_data->u)/(length);
	current_data->v = (v_t + current_data->v)/(length);
	current_data->x = (x_t + current_data->x)/(length);
	current_data->y = (y_t + current_data->y)/(length);
	current_data->score = MULTI*100.0*sqrt((current_data->u)*(current_data->u) + 
		(current_data->v)*(current_data->v))/sqrt((current_data->x)*
		(current_data->x) + (current_data->y)*(current_data->y));
	
	
	return current_data;
}


/*Searches through all array elements for array element with u flux closest to 
 *a percentage of the max u flux.
 *Requires data_t pointer to an array, an int size of array, double percentage
 *number between 0 and 1, and an output FILE pointer.
 *Returns the closes u flux value as a double
 */
double arraySearchFluxU(data_t* array[ARRAY_SIZE],int size
,double percentage,FILE *ofp){
	double max_flux = 0;
	double search_flux = 0;
	int i = 0;
//Finds max flux and percentage of max flux	
	 max_flux = array[size]->flux_u;
	 search_flux = percentage * max_flux;
	 
//Prints first search element	 
	 fprintf(ofp,"%.6lf",array[i]->flux_u);
	 i++;
	 
//Traverses array until it overshoots required flux	 
	 while(array[i]->flux_u <= search_flux){
	 	 fprintf(ofp,",%.6lf",array[i]->flux_u);
	 	 i++;
	 }

//Determines if overshot array entry is closer than previous undershot entry	 
	 if(fabs((array[i]->flux_u) - search_flux) < 
	 	 fabs((array[i-1]->flux_u) - search_flux)){
	 	fprintf(ofp,",%.6lf",array[i]->flux_u);
	 	return array[i]->flux_u;
	 }else{
	 	 return array[i-1]->flux_u;
	 }
	 return ERROR;
}


/*recursively searches through an array for a target value.
 *Requires data_t pointer to array,int lower search bound,int upper search 
 *bound, double search target value and an output FILE pointer.
 *Returns int TRUE if found, or int ERROR if not found
 */
int binarySearch(data_t* array[ARRAY_SIZE],int min, int max,
double target_flux_u,FILE *ofp){
	int mid = 0;
	
//Returns error if entry not found	
	if(min>max){
		return ERROR;
	}

//Finds midpoint between search bounds. If target less than mid, search again
//between min and mid and vice versa. Prints the midpoint to output file. 
	mid = min + (max - min) / 2;
	
	if(array[mid]->flux_u < target_flux_u){
		fprintf(ofp,"%.6lf,",array[mid]->flux_u);
		return binarySearch(array, mid+1, max, target_flux_u,ofp);
	}else if(array[mid]->flux_u > target_flux_u){
		fprintf(ofp,"%.6lf,",array[mid]->flux_u);
		return binarySearch(array, min, mid - 1, target_flux_u,ofp);
	}else{
		fprintf(ofp,"%.6lf",array[mid]->flux_u);
		return TRUE;
	}
	return ERROR;	
}


/*Searchest a linked list element wise for a target value.
 *Requires pointer to data_t, double target value and output FILE pointer as
 *arguments.
 *Returns int TRUE for success and int ERROR for failure.
 */
int searchLinkedList(data_t* data, double target_flux_u, FILE *ofp){
	
//Checks linked list until target value found. Prints all values to file	
	while(data->flux_u != target_flux_u){
		fprintf(ofp,"%.6lf,",data->flux_u);
		data = data->nextdata;
		
//If end of linked list is reached, return error
		if(!data){
			return ERROR;
		}
	}
	fprintf(ofp,"%.6lf",data->flux_u);
	return TRUE;
}

/*Searches through BST until target value is found.
 *Requires pointer to bst_t,double target value and output FILE pointer as
 *arguments.
 *Returns int TRUE if found and int ERROR if not found.
 */
int searchBST(bst_t* bst,double target_flux_u,FILE *ofp){
	node_t* current = bst->root;

//Checks current value against target value	and decides if it needs to 
//traverse left or right
	while(current->data->flux_u != target_flux_u){
		fprintf(ofp,"%.6lf,",current->data->flux_u);
		if(current->data->flux_u < target_flux_u){
			current = current->right;	
		}else if(current->data->flux_u > target_flux_u){
			current = current->left;	
		}
		if(!current){

//If end of branch is reached without findind the target value, returns ERROR
			return ERROR;
		}
	}
	fprintf(ofp,"%.6lf",current->data->flux_u);
	return TRUE;
}


/*Calculates w for a point(x,y).
 *Requires pointer to bst_t, v, and int array[] that has
 *the dimentions of the grid.
 */
void vortex(bst_t* bst,int x,int y,int xy[TWO_D]){
	double a = 0;
	double b = 0;
	double c = 0;
	double d = 0;
	data_t* data  = NULL;
	data_t* data_x  = NULL;
	data_t* data_y  = NULL;
	
//Calculates w depending on if the point falls into special case calculation
//puts w value into data_t struct->w
	if((x == (xy[X_COORD] - 1)) && ((y == xy[Y_COORD] - 1))){
		data = bstT(bst,x,y);
		data_x = bstT(bst,x-1,y); 
		data_y = bstT(bst,x,y-1);
		
		a = MULTI*(data->v) - (data_x->v);
		b = MULTI*(data->x) - (data_x->x);
		c = MULTI*(data->u) - (data_y->u);
		d = MULTI*(data->y) - (data_y->y);
		data->w = (a/b) - (c/d);
		
	}else if(x == (xy[X_COORD] - 1)){
		data = bstT(bst,x,y);
		data_x = bstT(bst,x-1,y); 
		data_y = bstT(bst,x,y+1);
		
		a = MULTI*(data->v) - (data_x->v);
		b = MULTI*(data->x) - (data_x->x);
		c = MULTI*(data_y->u) - (data->u);
		d = MULTI*(data_y->y) - (data->y);
		data->w = (a/b) - (c/d);
		
	}else if(y == (xy[Y_COORD] - 1)){
		data = bstT(bst,x,y);
		data_x = bstT(bst,x+1,y); 
		data_y = bstT(bst,x,y-1);
		
		a = MULTI*(data_x->v) - (data->v);
		b = MULTI*(data_x->x) - (data->x);
		c = MULTI*(data->u) - (data_y->u);
		d = MULTI*(data->y) - (data_y->y);
		data->w = (a/b) - (c/d);
	}else{
		data = bstT(bst,x,y);
		data_x = bstT(bst,x+1,y); 
		data_y = bstT(bst,x,y+1);
		
		a = MULTI*(data_x->v) - (data->v);
		b = MULTI*(data_x->x) - (data->x);
		c = MULTI*(data_y->u) - (data->u);
		d = MULTI*(data_y->y) - (data->y);
		data->w = (a/b) - (c/d);	
	}
	return;
}

/*Traverses bst to find data at a given x,y coord.
 *Requires pointer to bst_t, int x coord and in y coord as arguments.
 *Returns data_t pointer at that x,y location.
 */
data_t* bstT(bst_t* bst,int x,int y){
	node_t* current = bst->root;
	int i = 0;
	int j = 0;
	
//Traverses y branch for int y number of nodes	
	for(j = 0; j < y; j++){
		current = current->right;
	}
	
//Traverses x branch for int x number of nodes	
	for(i = 0; i < y; i++){
		current = current->left;
	}
	return current->data;
}


/*Prints data associated with a node_t struct.
 *Requires node_t pointer and output FILE pointer.
 */
void printNode(node_t* node,FILE *ofp){
	fprintf(ofp,"%.6lf,%.6lf,%.6lf,%.6lf,%.6lf\n",node->data->rho,node->data->u,
		node->data->v,node->data->x,node->data->y);
	return;
}

/*Prints data associated with a data_t struct.
 *Requires data_t pointer and output FILE pointer.
 */
void printDataS(data_t* data,FILE *ofp){
	fprintf(ofp,"%.6lf,%.6lf,%.6lf,%.6lf,%.6lf,%.6lf\n",data->rho
		,data->u,data->v,data->x,data->y,(data->score)/MULTI);
	return;
}


/*Gets a line from input FILE pointer and replaces the newline char with '\0'.
 *Requires char* string, int length of string and input FILE pointer.
 *Returns char* if success or fget's output if failure.
 */
char* fgetstr(char* string, int n, FILE *stream){
    char *result;
    result = fgets(string, n, stream);
    if(!result){
        return(result);
    }
//If newline char is found, replace it with string terminator    
    if(string[strlen(string) - 1] == '\n'){
        string[strlen(string) - 1] = '\0';
    }
    return(string);
}

/*Finds the dimetions of a 2D bst structure.
 *Requires pointer to bst_t and int array to put dimentions into.
 */
void findLength(bst_t* bst, int xy[TWO_D]){
	node_t* current = NULL;
	int x =0;
	int y = 0;
	
	current = bst->root;
	
//Finds length in y direction
	while(current->right != NULL){
		y++;
		current = current->right;
	}
	
//Finds length in x direction	
	while(current->left != NULL){
		x++;
		current = current->left;
	}
	xy[0] = x;
	xy[1] = y;
	return;
}

/*Gets the current time in microseconds.
 *Returns current time as a double.
 */
double getMicrotime(void){
	struct timeval currentTime;
	gettimeofday(&currentTime, NULL);
	return (double)currentTime.tv_sec * (double)1e6 + currentTime.tv_usec;
}

	
/*Frees all memory associated with a BST. 
 *Requires a bst_t pointer as an argument
 */
void freeBST(bst_t* bst,int data){
   freeTree(bst->root,data);
   free(bst);
}

/*Frees all memory associated with a BST recursively.
 *Requires a node_t pointer as an argument
 */
void freeTree(node_t* parent,int data){
    if(! parent){
    	return;
    }
    
//Calls freeTree for children branches then frees memory associated with 
//parent node.     
    freeTree(parent->left,data);
    freeTree(parent->right,data);
    if(data){
    freeData(parent->data);
    }
    free(parent);
    return;
}


/*Frees all memory associated with a data_t linked list recursively.
 *Requires a data_t pointer as an argument
 */
void freeData(data_t* parent){
	 if(! parent){
    	return;
    }
    
//Calls freeData for next data_t in list then frees memory associated with 
//parent data_t.     
    freeData(parent->nextdata);
    free(parent);
    return;
}
	
