function H = gen_hadamard_from_pairs(pairs)
M = length(pairs);
a = pairs{1};
N = length(a);
a = num2cell(a);

for i = 2:M
    aa = a;
    L = size(aa{1}, 1);
    a = cell(1, N);
    for j = 1:N
        a{j} = [aa{j}, zeros(L); zeros(L), aa{N-j+1}'];
    end

    bb = pairs{i};
    b = cell(1, N);
    for j = 1:N
        b{j} = [zeros(L), bb(j)*eye(L); -conj(bb(N-j+1))*eye(L), zeros(L)];
    end
    a = sp_sum(a, b);
end

CH = cell(N);
for i = 1:N
    CH(i, :) = [a(N-i+2:N), a(1:N-i+1)];
end
CH = cell2mat(CH);
H = zeros(2*size(CH));
for i = 1:size(CH, 1)
    for j = 1:size(CH, 2)
        H(2*i-1:2*i, 2*j-1:2*j) = sc2mat(CH(i, j));
    end
end

L = 2*size(a{1}, 1);
h = hadamard(L);
H = H*kron(eye(N), h);
end

function s = sp_sum(a, b)
N = length(a);
s = cell(1, N);
for i = 1:N
    s{i} = a{i}+b{i};
end
end

function mat = sc2mat(c)
if c==1
    mat = [1, 0; 0, 1];
elseif c==-1
    mat = [-1, 0; 0, -1];
elseif c==1i
    mat = [0, -1; 1, 0];
elseif c==-1i
    mat = [0, 1; -1, 0];
else
    mat = zeros(2);
end
end
