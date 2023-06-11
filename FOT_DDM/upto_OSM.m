mesh_spacing=0.05;
alpha=1/(mesh_spacing);
N=3;
load('partmat_OS')
partmat=OptimalSchwarz(partmat,alpha,N);