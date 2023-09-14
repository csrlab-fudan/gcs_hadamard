n1 = [3, 2]; % l_1, l_2, ...
n2 = [3, 4]; % m_1, m_2, ...
n = n1.*n2;
m = length(n);
seqs = cell(2*m, 1);
N = 4*sum(n);

old = 0;
for i = 1:m
    [~, ~, seqs1] = gen_quaternary_gcs_pair(n1(i));
    [~, ~, seqs2] = gen_quaternary_gcs_pair(n2(i));
    a = seqs1(1, :); b = seqs1(2, :);
    c = seqs2(1, :); d = seqs2(2, :);
    seqs{2*i-1} = [zeros(1, old), kron(a, c), zeros(1, N-2*old-2*n(i)), kron(b, d), zeros(1, old)];
    seqs{2*i} = [zeros(1, N/2-old-n(i)), kron(-conj(flip(a)), d), zeros(1, old*2), kron(conj(flip(b)), c), zeros(1, N/2-old-n(i))];
    old = old + n(i);
end

H = gen_hadamard_from_pairs(seqs);
result = H*H';