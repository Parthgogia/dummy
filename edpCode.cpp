#include <bits/stdc++.h>
using namespace std;
void solve(int row,int &n,vector<int>&arr, vector<vector<int>> &ans,
            vector<bool>&cols, vector<bool>&diag, 
            vector<bool>&antidiag){
            
            if(row==n) {ans.push_back(arr);return;}
            
            for (int i =0;i<n;i++){
                int ad = row+i;
                int d = row-i+n-1;
                
                if (cols[i] || diag[d] || antidiag[ad]) continue;
                    
                cols[i] = diag[d] = antidiag[ad] =true;
                arr.push_back(i+1);
                
                solve(row+1,n,arr,ans,cols,diag,antidiag);
                
                cols[i] = diag[d] = antidiag[ad] =false;
                arr.pop_back();
            } 
            
            }

int main() {
    int n =4;
    vector<bool> cols(n,false);
    vector<bool> diag(2*n-1,false);
    vector<bool> antidiag(2*n-1,false);
    vector<int> arr;
    vector<vector<int>> ans;

    solve(0,n,arr,ans,cols,diag,antidiag);

    for (auto it: ans){
        for (auto i: it){
            cout<<i<<" ";
        }
        cout<<endl;
    }

    return 0;
}





#include <bits/stdc++.h>
using namespace std;

bool solve(vector<vector<int>> &board,vector<vector<bool>> &cols,vector<vector<bool>> &rows,vector<vector<bool>> &boxes){

    for (int row = 0;row<9;row++){
        for (int col = 0;col<9;col++){
            if(board[row][col]==-1){
                int b = (row/3)*3 + col/3;
                for (int i =1;i<=9;i++){
                    if(!rows[row][i] && !cols[col][i] && !boxes[b][i]){

                        rows[row][i]=cols[col][i]=boxes[b][i]=true;
                        board[row][col] = i;

                        if(solve(board,cols,rows,boxes)) return true;

                        rows[row][i]=cols[col][i]=boxes[b][i]=false;
                        board[row][col] = -1;
                    }
                }
                return false;
            }
        }
    }
    return true;
}

int main() {
    vector<vector<bool>>rows(9,vector<bool>(10,false));
    vector<vector<bool>>cols(9,vector<bool>(10,false));
    vector<vector<bool>>boxes(9,vector<bool>(10,false));
    vector<vector<int>> board = {
    { 5,  3, -1,  -1,  7, -1,  -1, -1, -1},
    { 6, -1, -1,   1,  9,  5,  -1, -1, -1},
    {-1,  9,  8,  -1, -1, -1,  -1,  6, -1},

    { 8, -1, -1,  -1,  6, -1,  -1, -1,  3},
    { 4, -1, -1,   8, -1,  3,  -1, -1,  1},
    { 7, -1, -1,  -1,  2, -1,  -1, -1,  6},

    {-1,  6, -1,  -1, -1, -1,   2,  8, -1},
    {-1, -1, -1,   4,  1,  9,  -1, -1,  5},
    {-1, -1, -1,  -1,  8, -1,  -1,  7,  9}
    };
    
        for (int row = 0;row<9;row++){
        for (int col = 0;col<9;col++){
            if(board[row][col]!=-1){
                int num = board[row][col];
                int b = (row/3)*3 + col/3;
                rows[row][num] = cols[col][num] = boxes[b][num] = true;
            }
        }
    }

    solve(board,cols,rows,boxes);
    
    for (auto it:board){
        for (auto i:it){
            cout<<i<<" ";
        }
        cout<<endl;
    }


    return 0;
}


#include <bits/stdc++.h>
using namespace std;

// Function to check if it's safe to color the current vertex
// with the given color
bool issafe(int vertex, int col, vector<int> adj[], vector<int> &color) {
    for (auto it : adj[vertex]) {
        // If adjacent vertex has the same color, not safe
        if (color[it] != -1 && col == color[it])
            return false;
    }
    return true;
}

// Recursive function to try all colorings
bool cancolor(int vertex, int m, vector<int> adj[], vector<int> &color) {
    // If all vertices are colored successfully
    if (vertex == color.size())
        return true;

    // Try all colors from 0 to m-1
    for (int i = 0; i < m; i++) {
        if (issafe(vertex, i, adj, color)) {
            color[vertex] = i; 
            if (cancolor(vertex + 1, m, adj, color))
                // If the rest can be colored, return true
                return true; 
            color[vertex] = -1; 
        }
    }
    
    // No valid coloring found
    return false; 
}

bool graphColoring(int v, vector<vector<int>> &edges, int m) {
    vector<int> adj[v];

    // Build adjacency list from edges
    for (auto it : edges) {
        adj[it[0]].push_back(it[1]);
        adj[it[1]].push_back(it[0]); 
    }

    vector<int> color(v, -1); 
    return cancolor(0, m, adj, color);
}

int main() {
    int V = 4; 
    vector<vector<int>> edges = {{0, 1}, {0, 2},{0,3}, {1, 3}, {2, 3}}; 
    int m = 3; 

    // Check if the graph can be colored with m colors 
    // such that no adjacent nodes share the same color
    cout << (graphColoring(V, edges, m) ? "true" : "false") << endl;

    return 0;
}















// Online C++ compiler to run C++ program online
#include <bits/stdc++.h>
using namespace std;

void dfs(int node,vector<bool>&vis,stack<int>&st,vector<int> adj[]){
    vis[node] =true;
    for (auto it:adj[node]){
        if(!vis[it]) dfs(it,vis,st,adj);
    }
    st.push(node);
}

vector<int> toposort(vector<int> adj[],int n){
    stack<int>st;
    vector<bool> vis(n,0);
    for (int i =0;i<n;i++){
        if(!vis[i]) dfs(i,vis,st,adj);
    }
    vector<int> ans;
    while(!st.empty()){
        int node = st.top();st.pop();
        ans.push_back(node);
    }
    return ans;
}

int main(){
    vector<vector<int>> edges= {{2,3},{3,1},{4,0},{4,1},{5,0},{5,2}};
    int n =5;
    vector<int> adj[n+1];

    for (auto it:edges){
        adj[it[0]].push_back(it[1]);
    }

    // for(int i =0;i<=n;i++){
    //     cout<<i<<": ";
    //     for (auto it:adj[i]){
    //         cout<<it<<" ";
    //     }
    //     cout<<endl;
    // }
    auto sorted = toposort(adj,n+1);
    for (auto it:sorted){
        cout<<it<<" ";
    }
}











vector<int> kahnsort(vector<int> adj[],int n){
    queue<int>q;
    vector<int> indegree(n,0);
    for (int i =0;i<n;i++){
        for (auto it:adj[i]){
            indegree[it]++;
        }
    }

    for (int i=0;i<n;i++){
        if(indegree[i]==0) q.push(i);
    }
    vector<int> ans;
    while(!q.empty()){
        int node = q.front();q.pop();
        for (auto it:adj[node]){
            indegree[it]--;
            if(indegree[it]==0) q.push(it);
        }
        ans.push_back(node);
    }
    return ans;
}

int main(){
    vector<vector<int>> edges= {{2,3},{3,1},{4,0},{4,1},{5,0},{5,2}};
    int n =5;
    vector<int> adj[n+1];

    for (auto it:edges){
        adj[it[0]].push_back(it[1]);
    }

    // for(int i =0;i<=n;i++){
    //     cout<<i<<": ";
    //     for (auto it:adj[i]){
    //         cout<<it<<" ";
    //     }
    //     cout<<endl;
    // }
    auto sorted = kahnsort(adj,n+1);
    for (auto it:sorted){
        cout<<it<<" ";
    }
}



class Solution {
    private:
        void dfs(int node,stack<int>&st, vector<bool>& vis, vector<vector<int>> & adj){
            vis[node] =1;
            for (auto it:adj[node]){
                if(!vis[it]) dfs(it,st,vis,adj);
            }
            st.push(node);
        }
        void dfs2(int node,vector<bool>& vis, vector<vector<int>> & adj){
            vis[node] =1;
            for (auto it:adj[node]){
                if(!vis[it]) dfs2(it,vis,adj);
            }
        }
        
  public:
    int kosaraju(vector<vector<int>> &adj) {
        int n = adj.size();
        
        vector<bool> vis(n,0);
        stack<int> st;
        for (int i =0;i<n;i++){
            if(!vis[i]) dfs(i,st,vis,adj);
        }
        
        vector<vector<int>> adj2(n);
        for (int i =0;i<n;i++){
            for (auto it:adj[i]){
                adj2[it].push_back(i);
            }
        }       
        
        vis = vector<bool>(n,0);
        
        int scc=0;
        while(!st.empty()){
            int node = st.top();st.pop();
            if(!vis[node]){
                scc++;
                dfs2(node,vis,adj2);
            }
        }
        return scc;
        
    }
};



#include <bits/stdc++.h>
using namespace std;

// Utility to check if vertex v can be added at position pos in the Hamiltonian path
bool isSafe(int v, int pos,
            const vector<vector<int>>& graph,
            const vector<int>& path,
            const vector<bool>& used)
{
    // 1) Must be adjacent to the previous vertex in path
    if (!graph[ path[pos - 1] ][ v ]) 
        return false;

    // 2) Must not have been used already
    if (used[v]) 
        return false;

    return true;
}

bool hamCycleUtil(int pos,
                  const vector<vector<int>>& graph,
                  vector<int>& path,
                  vector<bool>& used)
{
    int V = graph.size();
    if (pos == V) {
        return graph[ path[pos - 1] ][ path[0] ] == 1;
    }

    for (int v = 1; v < V; v++) {
        if (isSafe(v, pos, graph, path, used)) {
            path[pos] = v;
            used[v] = true;

            if (hamCycleUtil(pos + 1, graph, path, used))
                return true;

            used[v] = false;
        }
    }
    return false;
}

bool hasHamiltonianCycle(const vector<vector<int>>& graph) {
    int V = graph.size();
    vector<int> path(V, -1);
    vector<bool> used(V, false);

    path[0] = 0;
    used[0] = true;

    if (hamCycleUtil(1, graph, path, used)) {
        cout << "Hamiltonian cycle found: ";
        for (int v : path) 
            cout << v << " ";
        cout << path[0] << "\n"; // return to start
        return true;
    }

    cout << "No Hamiltonian cycle exists\n";
    return false;
}

int main() {
    // Example graph (adjacency matrix)
    //   (0)--(1)
    //    |    |
    //   (3)--(2)
    vector<vector<int>> graph = {
        {0, 1, 0, 1},
        {1, 0, 1, 0},
        {0, 1, 0, 1},
        {1, 0, 1, 0}
    };

    hasHamiltonianCycle(graph);
    return 0;
}

#include <bits/stdc++.h>
using namespace std;

// Search the pat string in the txt string 
void search(string pat, string txt, int q)
{
    int M = pat.size();
    int N = txt.size();
    int i, j;
    int p = 0; // hash value for pattern
    int t = 0; // hash value for txt
    int h = 1;
    int d = 256; // d is the number of characters in the input alphabet 
    // The value of h would be "pow(d, M-1)%q"
    for (i = 0; i < M - 1; i++)
        h = (h * d) % q;

    // Calculate the hash value of pattern and first
    // window of text
    for (i = 0; i < M; i++) {
        p = (d * p + pat[i]) % q;
        t = (d * t + txt[i]) % q;
    }

    // Slide the pattern over text one by one
    for (i = 0; i <= N - M; i++) {

        // Check the hash values of current window of text
        // and pattern. If the hash values match then only
        // check for characters one by one
        if (p == t) {
            /* Check for characters one by one */
            for (j = 0; j < M; j++) {
                if (txt[i + j] != pat[j]) {
                    break;
                }
            }

            // if p == t and pat[0...M-1] = txt[i, i+1,
            // ...i+M-1]

            if (j == M)
                cout << "Pattern found at index " << i
                     << endl;
        }

        // Calculate hash value for next window of text:
        // Remove leading digit, add trailing digit
        if (i < N - M) {
            t = (d * (t - txt[i] * h) + txt[i + M]) % q;

            // We might get negative value of t, converting
            // it to positive
            if (t < 0)
                t = (t + q);
        }
    }
}

/* Driver code */
int main()
{
    string txt = "GEEKS FOR GEEKS";
    string pat = "GEEK";

    // we mod to avoid overflowing of value but we should
    // take as big q as possible to avoid the collison
    int q = INT_MAX;
    // Function Call
    search(pat, txt, q);
    return 0;
}





#include<bits/stdc++.h>
 
using namespace std;
int kmp(string String, string pattern) {
  int i = 0, j = 0, m = pattern.length(), n = String.length();
  pattern = ' ' + pattern; //just shifting the pattern indices by 1
  vector < int > piTable(m + 1, 0);
  for (int i = 2; i <= m; i++) {
    while (j <= m && pattern[j + 1] == pattern[i])
      piTable[i++] = ++j;
    j = 0;
  }
  j = 0;
  for (int i = 0; i < n; i++) {
    if (pattern[j + 1] != String[i]) {
      while (j != 0 && pattern[j + 1] != String[i])
        j = piTable[j];
    }
    j++;
    if (j == m) return i - m + 1;
  }
  return -1;

}
int main() {
  string pattern="aaaaaab", String="aaaaaaaamaaaaaab";

  int index = kmp(String, pattern);
  if (index == -1) cout << "The pattern is not found";
  else cout << "The pattern " << pattern << " is found in the given string " 
  << String << " at " << index;
  return 0;
}

