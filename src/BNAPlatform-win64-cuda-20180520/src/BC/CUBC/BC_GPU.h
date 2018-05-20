#define NODE 1
#define EDGE 0

#define ARRAY 0
#define PRED 1
#define SUCC 0

void Betweenness_GPU_node_array(int *r, int *c, int numVertices, int numEdges, float *BC, int grid, int thread);
void Betweenness_GPU_node_pred(int *r, int *c, int numVertices, int numEdges, float *BC, int grid, int thread);
void Betweenness_GPU_node_succ(int *r, int *c, int numVertices, int numEdges, float *BC, int grid, int thread);
void Betweenness_GPU_edge_array(int *r, int *r_full, int *c, int numVertices, int numEdges, float *BC, int grid, int thread);
void Betweenness_GPU_edge_pred(int *r, int *r_full, int *c, int numVertices, int numEdges, float *BC, int grid, int thread);
void Betweenness_GPU_edge_succ(int *r, int *r_full, int *c, int numVertices, int numEdges, float *BC, int grid, int thread);


