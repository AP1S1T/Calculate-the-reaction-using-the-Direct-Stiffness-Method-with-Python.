import numpy as np

def input_list(function):
    return list(map(int,input(function).split()))

n = int(input("Input number of members:")) #number of members
I = input_list("Enter the moment of inertia for each member:") #moment of Inertia
L = input_list("Enter the length for each member:") #length in m
Du = int(input("Enter the number of unconstained degrees of freedom:"))
Dc = int(input("Enter the number of constained degrees of freedom:"))
dof = Du+Dc  #total degrees of freedom
Dul = input_list("Enter the indices of Unconstained degrees of freedom :")
Dcl = input_list("Enter the indices of constained degrees of freedom :")

l = [] #empty list for all global labels 

for i in range(n):
    member_label = f"member{i+1}"
    l.extend(input_list(f"Enter the local labels for {member_label}:"))

Ktotal = np.zeros((dof,dof)) #Globals stiffness matrix size (6x6)
for i in range(n):
    member_label = f"member{i+1}"
    fem = input_list(f"Enter the local fix-end moment for{member_label}:")

# rotation coefficients for each member
    rc1 = 12 * np.array(I)/np.array(L)**3
    rc2 = 6 * np.array(I)/np.array(L)**2
    rc3 = 4 * np.array(I)/np.array(L)
    rc4 = 2 * np.array(I)/np.array(L)

 #Stiffness matrix 4x4 
for i in range(n):
    knew = np.zeros((dof, dof))
    k1 = np.array([rc1[i], rc2[i], -rc1[i], rc2[i]])
    k2 = np.array([rc2[i], rc3[i], -rc2[i], rc4[i]])
    k3 = np.array([-rc1[i], -rc2[i], rc1[i], -rc2[i]])
    k4 = np.array([rc2[i], rc4[i], -rc2[i], rc3[i]])
    K = np.array([k1, k2, k3, k4])
    print(f'Member Numbers = {i + 1}')
    print('Local stiffness matrix of member, [K] = ')
    print(K)

    local_indices = l[i * 4: (i + 1) * 4]

    print(f'Local indices for member {i + 1}: {local_indices}')

    for p in range(4):
        for q in range(4):
            global_row = local_indices[p]-1
            global_col = local_indices[q]-1
            knew[global_row,global_col] += K[p,q]
    
    Ktotal += knew
print('Stiffness Matrix of complete structure, [Ktotal] = ')
print(Ktotal)

Kunc = np.zeros((Du,Du)) #Unconstrained matrix

for x in range(Du):
    for y in range(Du):
        Kunc[x,y] = Ktotal[x,y]

print('Unconstrained Stiffness Matrix, [Kunc] = ')
print(Kunc)

KuncInv = np.linalg.inv(Kunc)
print('Inverse Unconstrained stiffness matrix = ')
print(KuncInv)

Qk = np.array(input_list('Input your Nodal Loads ='))
delu = np.dot(KuncInv,Qk)
print('Displacement/EI=',delu)
dec = np.zeros((Dc,1))

Kcon = np.zeros((Dc,Du))

for x in range(Dc):
    for y in range(Du):
        Kcon[x,y] = Ktotal[Du+x,y]
print('Constrained stiffness matrix,[Kcon] = ')
print(Kcon)

Qu = np.dot(Kcon,delu)
for i, value in enumerate(np.round(Qu.flatten(), 3), start=Du+1):
    print(f"Q{i} =", value)

