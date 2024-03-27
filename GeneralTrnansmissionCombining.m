function GeneralTrnansmissionCombining(n_t, n_r)
rand('seed',1);
n_t=2;n_r=3;
plotDebugFlag=1;

n=n_t*n_r;
% n=9;
min_p = find_minimum_p(n);
hr_set = build_hurwitz_radon(min_p);
mat_set = hr_set(:,:,1:n-1);
mat_set(:,:,size(mat_set,3)+1)=eye(min_p);
% Construct variables vector
% var_vector = zeros(min_p, 1);
% basic_sym_string = 'x_';
% result_matrix = zeros(min_p, n);
% for i in range(p):
% for i=1:p
%         %sym = Symbol(basic_sym_string + str(i + 1))
%         sym{i}=[sym basic_sym_string+num2str(i)];
%         %var_vector[i, 0] = sym
% end
xVectReal=sym('x%d', [min_p 1]);
% xVectComplex=xVectReal;
% for kk=1:length(xVectComplex)
%     xVectComplex(kk)=conj(xVectComplex(kk));
% end
hVect=sym('h%d', [n 1]);
% Construct OSTBC
%for i in range(n):
result_matrixReal=sym('y%d', [min_p n]);
% result_matrixComplex=sym('y%d', [min_p n]);
for i=1:n
    % result_matrix[:, i] = mat_set[i] * var_vector;
    [rowInd,colInd]=find(mat_set(:,:,i)~=0);
    if plotDebugFlag
        disp(i);
    end
    for kk=1:min_p        
        result_matrixReal(rowInd(kk),end-i+1)=mat_set(rowInd(kk),colInd(kk),i)*xVectReal(colInd(kk));
    end

    % result_matrixReal(:, end-i+1) = mat_set(:,:,i) * xVectReal;
    % result_matrixComplex(:, end-i+1) = mat_set(:,:,i) * xVectComplex;
end
%OSTBCmatrix=[result_matrixReal;result_matrixComplex];
OSTBCmatrix=[result_matrixReal]
OSTBCafterChannel=OSTBCmatrix*hVect;
if plotDebugFlag
    OSTBCafterChannel
end
transmitMatrixMIMO=sym(zeros(n_t*min_p,n_t));
for kk=1:n_t
    transmitMatrixMIMO((kk-1)*min_p+(1:min_p),kk)=xVectReal;
end
hMat=reshape(hVect,n_t,n_r);
yReal=transmitMatrixMIMO*hMat;
if plotDebugFlag
    yReal
end


for rowInd=1:size(OSTBCmatrix,1)
    if plotDebugFlag
        disp(rowInd);
    end
    % opernads=children(OSTBCafterChannel(rowInd));
    for colInd=1:size(OSTBCmatrix,2)
        value=OSTBCmatrix(rowInd,colInd)*hVect(colInd);
        % [col,row]=find(opernads{colInd}==yReal);
        [col,row]=find(value==yReal,1);
        if ~isempty(col) 
            if colInd==1
                combMat(rowInd,1)=sym(['y_' num2str(col) '_' num2str(row)]);
            else
                combMat(rowInd,1)=combMat(rowInd,1)+sym(['y_' num2str(col) '_' num2str(row)]);
            end
        end
        % [col,row]=find(opernads{colInd}==-yReal);
        [col,row]=find(value==-yReal,1);
        if ~isempty(col) 
            if colInd==1
                combMat(rowInd,1)=-sym(['y_' num2str(col) '_' num2str(row)]);
            else
                combMat(rowInd,1)=combMat(rowInd,1)-sym(['y_' num2str(col) '_' num2str(row)]);
            end
        end
    end
end
transmitMatrixMIMO
combMat
keyboard;
end

function calc_result=calc_number_of_iters(n_t, n_r)
%"""Calculates mathematically the number of iterations of the generator 'generate_all_matrix_options"""
calc_result = 1;
for row_ind=0:(n_r - 2)
    calc_result = calc_result * nchoosek(n_t * (n_r - row_ind), n_t) * factorial(n_t) * (2 ^ (n_t-1));
end
end

function p=find_minimum_p(n)
%"""The function find for 'n' number of transmit antennas the minimum delay 'p'"""
p = 1;
while (1)
    if (hurwitz_radon_number(p) >= n)
        break;
    end
    p=p+1;
end
end

function out=hurwitz_radon_number(n)
%""" The functions gets as input a number 'n' and returns its Hurwitz-Radon Number raw(n)"""
[a, b, c, d] = find_hurwitz_radon_parameters(n);
out=8*c + 2^d;
end

function [a, b, c, d]=find_hurwitz_radon_parameters(n)
%""" The functions gets as input a number 'n', and returns the parameters (a,b,c,d) that associated
%with the number. Those numbers are needed for Hurwitz-radon number computation."""
for a=0:n-1
    if mod(n,(2^a)) ~= 0
        continue;
    else
        if mod(idivide(int64(n),int64(2 ^ a),'floor'),2) == 1  %# is odd
            b = idivide(int64(n),int64(2 ^ a),'floor');
            break;
        end
    end
end
c = idivide(int64(a),int64(4),'floor');
d=mod(a,4);
end

function mat_set=build_hurwitz_radon(n)
%""" The function gets as input a number 'n' and returns a set of matrices of size n x n
%that meets the Hurwitz-Radon constraints. The functions assumes that 'b=1', where b is a Hurwitz-Radon parameter."""
mat_set=[];
if (n == 0)
    return;
else
    [a, b, c, d] = find_hurwitz_radon_parameters(n);
    s = c - 1;
    R = [[0 1];[-1 0]];%Matrix([[0, 1], [-1, 0]])
    P = [[0 1];[1 0]];%Matrix([[0, 1], [1, 0]])
    Q = [[1 0];[0 -1]];%Matrix([[1, 0], [0, -1]])
    I = eye(2);
    %mat_set = [];
    if (a == 0)
        return;
    elseif (a == 1)
        mat_set=R;
%         if size(mat_set,1)==0
%             mat_set=R;
%         else
%             mat_set(:,:,size(mat_set,3)+1)=R;
%         end
        %mat_set=[mat_set R];
        return;
    elseif (a == 2)
        mat_set(:,:,1)=kron(R, I);
        mat_set(:,:,2)=kron(P, R);
        mat_set(:,:,3)=kron(Q, R);
        %mat_set = [kron(R, I), kron(P, R), kron(Q, R)];
        return;
    elseif (a == 3)
        mat_set(:,:,1)=kron(I, kron(R, I));
        mat_set(:,:,2)=kron(I, kron(P, R));
        mat_set(:,:,3)=kron(Q, kron(Q, R));
        mat_set(:,:,4)=kron(P, kron(Q, R));
        mat_set(:,:,5)=kron(R, kron(P, Q));
        mat_set(:,:,6)=kron(R, kron(P, P));
        mat_set(:,:,7)=kron(R, kron(Q, I));
        %         mat_set = [kron(I, kron(R, I)),
        %             kron(I, kron(P, R)),
        %             kron(Q, kron(Q, R)),
        %             kron(P, kron(Q, R)),
        %             kron(R, kron(P, Q)),
        %             kron(R, kron(P, P)),
        %             kron(R, kron(Q, I))];
        return;
    elseif(d == 0)
        %mat_set.append(kron(R, eye(2 ** (4*s + 3))))
        %             mat_set=[mat_set kron(R, eye(2 ^ (4*s + 3)))];
        if size(mat_set,1)==0
            mat_set=kron(R, eye(2 ^ (4*s + 3)));
        else
            mat_set(:,:,size(mat_set,3)+1)=kron(R, eye(2 ^ (4*s + 3)));
        end
        %set_A = build_hurwitz_radon(2 ** (4*s + 3))
        set_A = build_hurwitz_radon(2 ^(4*s + 3));
        for matInd=1:size(set_A,3)
            %                 mat_set=[mat_set kron(Q,set_A(:,:,matInd))];
            mat_set(:,:,size(mat_set,3)+1)=kron(Q,set_A(:,:,matInd));
        end
        %             for A in set_A:
        %                  mat_set.append(kron(Q, A))
        return;
    else  %# if 1<= d < 4
        set_L = build_hurwitz_radon(2 ^ (4*s + 3));
        set_A = build_hurwitz_radon(2 ^ d);
        for matInd=1:size(set_A,3)
            %                 mat_set=[mat_set kron(P,kron(eye(2^(4*s+3)),set_A(matInd,:,:)))];
            if size(mat_set,1)==0
                mat_set=kron(P,kron(eye(2^(4*s+3)),set_A(:,:,matInd)));
            else
                mat_set(:,:,size(mat_set,3)+1)=kron(P,kron(eye(2^(4*s+3)),set_A(:,:,matInd)));
            end
        end
        %             for A in set_A:
        %                 mat_set.append(kron(P, kron(eye(2 ** (4*s + 3)), A)))
        for matInd=1:size(set_L,3)
            %                 mat_set=[mat_set kron(Q,kron(set_L(matInd,:,:)),eye(2^d))];
            mat_set(:,:,size(mat_set,3)+1)=kron(Q,kron(set_L(:,:,matInd),eye(2^d)));
            %             for L in set_L:
            %                 mat_set.append(kron(Q, kron(L, eye(2 ** d))))
        end
        mat_set(:,:,size(mat_set,3)+1)=kron(R, eye(2 ^ (4*s+3) * 2^d));
        return;
    end
end
end

function sigMat=generateSigMat(n_t, n_r)
numOfindices=n_t*n_r-(n_t+n_r-1);
combinations=dec2bin(0:2^numOfindices-1);
relIndices=[];
for ind=1:n_t*n_r
    if mod(ind,n_t)~=1 && mod(ind,n_r)~=0
        relIndices=[relIndices ind];
    end
end
for matInd=1:length(combinations)
    tempMat=zeros(n_r,n_t);
    for charInd=1:length(relIndices)
        tempMat(relIndices(charInd))=str2num(combinations(matInd,charInd));
    end
    sigMat(:,:,matInd)=tempMat;
end
end

function permMat=generatePermMat(n_t, n_r)
numOfindices=n_t*n_r;
if factorial(numOfindices)<1e7
    combinations=perms(1:numOfindices);
    for permInd=1:length(combinations)
        permMat(:,:,permInd)=reshape(combinations(permInd,:),n_r,n_t);
    end
else
    permMat=zeros(n_r,n_t,1e7);
    for permInd=1:1e7
        permMat(:,:,permInd)=reshape(randperm(numOfindices),n_r,n_t);
    end
end
end


% function generate_block_matrix(mat_set, n_t, n_r)
% for mat_ind, sign_mat in generate_all_matrix_options(n_t, n_r):
%     big_mat = Matrix([[mat_set[ind] * (sign_mat[row][col] * 2 - 1) for col, ind in enumerate(indices)]
%         for row, indices in enumerate(mat_ind)])
%             yield big_mat
%
%             function generate_all_matrix_options(n_t, n_r)
%                 %""" generates all possible matrices out of the subset with one exception:
%                 %the 'last_row' row contains the only option for values without checking their permutations.
%                 %For each such matrix it also generates all possible sign patterns ('0' or '1') whereas the first column
%                 %and the last row are fixed to '0'."""
%                 %# run over options from first row to last
%                 for x, y in generate_from_row(0, list(range(n_t * n_r)), n_t, n_r - 1):
%                     yield x, y

% function [sub_result, sub_pattern]=generate_from_row(row_num, subset, n_t, last_row)
% %""" generates all possible matrices out of the subset, i.e. going through all over possible
% %options for row number 'row_num' and all possible options for the rows below with one exception:
% %the 'last_row' row contains the only option for values without checking their permutations.
% %For each such matrix it also generates all possible sign patterns ('0' or '1') whereas the first column
% %and the last row are fixed to '0'."""
% if row_num == last_row
%     %sub_pattern=[[0]*n_t];
%     sub_pattern=zeros(1,3);
%     return;
% end
% for in_, out_ in comb_and_comp(subset, n_t):
%         # run over all possible permutations
%         for permut in itertools.permutations(in_, n_t):
%             for sign_pattern in generate_sign_pattern(n_t):
%                 g = generate_from_row(row_num + 1, out_, n_t, last_row)
%                 for sub_result, sub_pattern in g:
%                     sub_result.insert(0, list(permut))
%                     sub_pattern.insert(0, list(sign_pattern))
%                     yield sub_result, sub_pattern

% function [in_,out_]=comb_and_comp(lst, n)
% %     %""" This code was taken from 'Stack Overflow'
% %     %The generator returns all possible sub-list of size 'n' and their complements as  a tuple of lists.
% %     %i.e returns: (sub_list, complement)"""
% %     %# no combinations
% %     if (len(lst)) < n
% %         return
% %     end
% %     %# trivial 'empty' combination
% %     if (n == 0 || lst == [])
% % %         yield [], lst
% %     else
% % %         first, rest = lst[0], lst[1:]
% %         first = lst[0];
% %         rest = lst[1:end];
% %         %# combinations that contain the first element
% %         [in_temp,out_temp]=comb_and_comp(lst, n);
% %         for ind=1:length()
% %         for in_, out in comb_and_comp(rest, n - 1):
% %             yield [first] + in_, out
% %         # combinations that do not contain the first element
% %         for in_, out in comb_and_comp(rest, n):
% %             yield in_, [first] + out
% %             end
% all_perms=prems(1:length(lst));
% in_= lst(all_perms(1:n));
% out_=lst(all_perms(n+1:end));
% end
