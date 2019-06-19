function [ ] = bchverilog(k, t, make_testbench)
%bchverilog Creates Verilog codec for a shortened systematic BCH code

% Find the code parameters
m = 2;
T = zeros(1, 3);
a = 1;
while T(a, 2) < k
    m = m+1;
    T = bchnumerr(2^m-1);
    a = 1;
    while a < size(T, 1) && T(a, 3) < t
        a = a+1;
    end
end
l = T(a, 2)-k;
n = T(a, 1)-l;
fid = fopen('n.txt', 'w');
fprintf(fid, '%u\n', n);
fclose(fid);

% Write the Verilog for the encoder
genpoly = bchgenpoly(n+l, k+l);
[~, G] = cyclgen(n+l, fliplr(double(genpoly.x)), 'system');
G = G(1:k, 1:n);
fid = fopen('bch_encoder.v', 'w');
fprintf(fid, 'module bch_encoder(message, codeword);\n');
fprintf(fid, '  input [%u:0] message;\n', k-1);
fprintf(fid, '  output [%u:0] codeword;\n\n', n-1);
for a = 1:n
    fprintf(fid, '  assign codeword[%u] = ', a-1);
    if ~sum(G(:, a))
        fprintf(fid, '1''b0\n');
    else
        c = 0;
        for b = 1:k
            if G(b, a)
                if c
                    fprintf(fid, ' ^ ');
                end
                c = c+1;
                fprintf(fid, 'message[%u]', b-1);
            end
        end
        fprintf(fid, ';\n');
    end
end
fprintf(fid, 'endmodule\n');
fclose(fid);

% Write the Verilog for the syndromes
H = zeros(2*t-1, n, m);
for a = 1:2*t-1
    H(a, :, :) = fliplr(gftuple(a*(0:n-1).', m));
end
fid = fopen('bch_syndromes.v', 'w');
fprintf(fid, 'module bch_syndromes(received');
for a = 1:2*t-1
    fprintf(fid, ', syndrome%u', a);
end
fprintf(fid, ');\n');
fprintf(fid, '  input [%u:0] received;\n', n-1);
for a = 1:2*t-1
    fprintf(fid, '  output [%u:0] syndrome%u;\n', m-1, a);
end
fprintf(fid, '\n');
for a = 1:2*t-1
    fprintf(fid, '  assign syndrome%u = ', a);
    for b = 0:n-1
        if b
            fprintf(fid, ' ^ ');
        end
        fprintf(fid, '({%u{received[%u]}} & %u''b%s)', m, b, m, ...
            num2str(H(a, b+1, :), '%u'));
    end
    fprintf(fid, ';\n');
end
fprintf(fid, 'endmodule\n');
fclose(fid);

% Write the Verilog for the Galois field multiplier
fid = fopen('bch_gfmult.v', 'w');
fprintf(fid, 'module gfmult(a, b, c);\n');
fprintf(fid, '  input [%u:0] a;\n', m-1);
fprintf(fid, '  input [%u:0] b;\n', m-1);
fprintf(fid, '  output [%u:0] c;\n\n', m-1);
fprintf(fid, '  assign c = ');
for a = 0:m-1
    if a
        fprintf(fid, ' ^ ');
    end
    fprintf(fid, '({%u{b[%u]}} & ((a << %u)', m, a, a);
    if a
        T = fliplr(gftuple(m+a-(1:a).', m));
    else
        T = fliplr(gftuple(m+a-(1:-1:0).', m));
    end
    for b = 1:a
        fprintf(fid, ' ^ ({%u{a[%u]}} & %u''b%s)', m, m-b, m, ...
            num2str(T(b,:), '%u'));
    end
    fprintf(fid, '))');
end
fprintf(fid, ';\n');
fprintf(fid, 'endmodule\n');
fclose(fid);

% Write the Verilog for the inversionless Berlekamp-Massey algorithm
gfzero = fliplr(gftuple(0, m));
fid = fopen('bch_ibm.v', 'w');
fprintf(fid, 'module bch_ibm(');
for a = 1:2*t-1
    fprintf(fid, 'syndrome%u, ', a);
end
for a = 0:t
    if a
        fprintf(fid, ', ');
    end
    fprintf(fid, 'locator%u', a);
end
fprintf(fid, ');\n');
for a = 1:2*t-1
    fprintf(fid, '  input [%u:0] syndrome%u;\n', m-1, a);
end
for a = 0:t
    fprintf(fid, '  output [%u:0] locator%u;\n', m-1, a);
end
fprintf(fid, '\n');
% Declare the wires
for a = 0:t
    for b = 0:t
        fprintf(fid, '  wire [%u:0] nu%u_%u;\n', m-1, 2*a, b);
    end
end
for a = 0:t
    for b = 0:t
        fprintf(fid, '  wire [%u:0] kappa%u_%u;\n', m-1, 2*a, b);
    end
end
for a = 0:t-1
    fprintf(fid, '  wire [%u:0] delta%u;\n', m-1, 2*a);
end
for a = 0:t-1
    fprintf(fid, '  wire [%u:0] d%u;\n', m-1, 2*a);
end
for a = 1:t-1
    for b = 0:t
        fprintf(fid, '  wire [%u:0] delta%u_x_nu%u_%u;\n', ...
            m-1, 2*(a-1), 2*a, b);
    end
end
for a = 0:t-1
    for b = 0:t
        fprintf(fid, '  wire [%u:0] d%u_x_kappa%u_%u;\n', ...
            m-1, 2*a, 2*a, b);
    end
end
for a = 0:t-1
    c = 2*a+1;
    if c > t
        c = t;
    end
    for b = 0:c
        if 2*a+1-b > 0
            fprintf(fid, '  wire [%u:0] nu%u_%u_x_syndrome%u;\n', ...
                m-1, 2*a, b, 2*a+1-b);
        end
    end
end
for a = 0:t-1
    fprintf(fid, '  wire cond%u;\n', 2*a);
end
fprintf(fid, '\n');
% Instantiate the Galois field multipliers
for a = 1:t-1
    for b = 0:t
        fprintf(fid, '  gfmult gfm_delta%u_x_nu%u_%u(', ...
            2*(a-1), 2*a, b);
        fprintf(fid, 'delta%u, nu%u_%u, ', 2*(a-1), 2*a, b);
        fprintf(fid, 'delta%u_x_nu%u_%u);\n', 2*(a-1), 2*a, b);
    end
end
for a = 0:t-1
    for b = 0:t
        fprintf(fid, '  gfmult gfm_d%u_x_kappa%u_%u(', 2*a, 2*a, b);
        fprintf(fid, 'd%u, kappa%u_%u, ', 2*a, 2*a, b);
        fprintf(fid, 'd%u_x_kappa%u_%u);\n', 2*a, 2*a, b);
    end
end
for a = 0:t-1
    c = 2*a+1;
    if c > t
        c = t;
    end
    for b = 0:c
        if 2*a+1-b > 0
            fprintf(fid, '  gfmult gfm_nu%u_%u_x_syndrome%u(', ...
                2*a, b, 2*a+1-b);
            fprintf(fid, 'nu%u_%u, syndrome%u, ', 2*a, b, 2*a+1-b);
            fprintf(fid, 'nu%u_%u_x_syndrome%u);\n', 2*a, b, 2*a+1-b);
        end
    end
end
fprintf(fid, '\n');
% Assign the wires
for a = 0:t
    fprintf(fid, '  assign locator%u = nu%u_%u;\n', a, 2*t, a);
end
fprintf(fid, '  assign nu0_0 = %u''b%s;\n', m, ...
    num2str(gfzero, '%u'));
for a = 1:t
    fprintf(fid, '  assign nu0_%u = %u''b0;\n', a, m);
end
for a = 1:t
    for b = 0:t
        if a > 1
            fprintf(fid, '  assign nu%u_%u = delta%u_x_nu%u_%u', ...
                2*a, b, 2*a-4, 2*a-2, b);
        else
            fprintf(fid, '  assign nu%u_%u = nu%u_%u', ...
                2*a, b, 2*a-2, b);
        end
        if b
            fprintf(fid, ' ^ d%u_x_kappa%u_%u', 2*a-2, 2*a-2, b-1);
        end
        fprintf(fid, ';\n');
    end
end
fprintf(fid, '  assign kappa0_0 = %u''b%s;\n', m, ...
    num2str(gfzero, '%u'));
for a = 1:t
    fprintf(fid, '  assign kappa0_%u = %u''b0;\n', a, m);
end
for a = 1:t
    for b = 0:t
        fprintf(fid, '  assign kappa%u_%u = cond%u ? ', 2*a, b, 2*a-2);
        if b < 2
            fprintf(fid, '%u''b0 : ', m);
        else
            fprintf(fid, 'kappa%u_%u : ', 2*a-2, b-2);
        end
        if b < 1
            fprintf(fid, '%u''b0;\n', m);
        else
            fprintf(fid, 'nu%u_%u;\n', 2*a-2, b-1);
        end
    end
end
for a = 0:t-1
    fprintf(fid, '  assign delta%u = cond%u ? ', 2*a, 2*a);
    if ~a
        fprintf(fid, '%u''b%s', m, num2str(gfzero, '%u'));
    else
        fprintf(fid, 'delta%u', 2*a-2);
    end
    fprintf(fid, ' : d%u;\n', 2*a);
end
for a = 0:t-1
    fprintf(fid, '  assign d%u = ', 2*a);
    c = 2*a+1;
    if c > t
        c = t;
    end
    for b = 0:c
        if b
            fprintf(fid, ' ^ ');
        end
        if 2*a+1-b > 0
            fprintf(fid, 'nu%u_%u_x_syndrome%u', 2*a, b, 2*a+1-b);
        else
            fprintf(fid, 'nu%u_%u', 2*a, b);
        end
    end
    fprintf(fid, ';\n');
end
for a = 0:t-1
    fprintf(fid, '  assign cond%u = (~| d%u)', 2*a, 2*a);
    for b = a+1:t
        fprintf(fid, ' | (| nu%u_%u)', 2*a, b);
    end
    fprintf(fid, ';\n');
end
fprintf(fid, 'endmodule\n');
fclose(fid);

% Write the Verilog for the Chien search and correction
fid = fopen('bch_chien.v', 'w');
fprintf(fid, 'module bch_chien(received');
for a = 0:t
    fprintf(fid, ', locator%u', a);
end
fprintf(fid, ', codeword);\n');
fprintf(fid, '  input [%u:0] received;\n', n-1);
for a = 0:t
    fprintf(fid, '  input [%u:0] locator%u;\n', m-1, a);
end
fprintf(fid, '  output [%u:0] codeword;\n\n', n-1);
fprintf(fid, '  wire [%u:0] error;\n', n-1);
for a = 0:n-1
    for b = 0:t
        fprintf(fid, ...
            '  wire [%u:0] alpha%u_%u_x_locator%u;\n', m-1, a, b, b);
    end
end
fprintf(fid, '\n');
for a = 0:n-1
    T = fliplr(gftuple((n+l-a)*(0:t).', m));
    for b = 0:t
        fprintf(fid, '  gfmult gfm_alpha%u_%u_x_locator%u(', a, b, b);
        fprintf(fid, '%u''b%s, ', ...
            m, num2str(T(b+1,:), '%u'));
        fprintf(fid, 'locator%u, alpha%u_%u_x_locator%u);\n', b, a, b, b);
    end
end
fprintf(fid, '\n');
fprintf(fid, '  assign codeword = received ^ error;\n');
fprintf(fid, '  assign error = {');
for a = n-1:-1:0
    if a < n-1
        fprintf(fid, ', ');
    end
    fprintf(fid, '(~| (');
    for b = 0:t
        if b
            fprintf(fid, ' ^ ');
        end
        fprintf(fid, 'alpha%u_%u_x_locator%u', a, b, b);
    end
    fprintf(fid, '))');
end
fprintf(fid, '};\n');
fprintf(fid, 'endmodule\n');
fclose(fid);

% Write the Verilog for the detector
fid = fopen('bch_detector.v', 'w');
fprintf(fid, 'module bch_detector(received, ');
for a = 1:2*t-1
    fprintf(fid, 'syndrome%u, ', a);
end
fprintf(fid, 'detection);\n');
fprintf(fid, '  input [%u:0] received;\n', n-1);
for a = 1:2*t-1
    fprintf(fid, '  output [%u:0] syndrome%u;\n', m-1, a);
end
fprintf(fid, '  output detection;\n\n');
fprintf(fid, '  bch_syndromes syndromes(received');
for a = 1:2*t-1
    fprintf(fid, ', syndrome%u', a);
end
fprintf(fid, ');\n');
fprintf(fid, '  assign detection = ');
for a = 1:2*t-1
    if a > 1
        fprintf(fid, ' | ');
    end
    fprintf(fid, '(| syndrome%u)', a);
end
fprintf(fid, ';\n');
fprintf(fid, 'endmodule\n');
fclose(fid);

% Write the Verilog for the corrector
fid = fopen('bch_corrector.v', 'w');
fprintf(fid, 'module bch_corrector(received, ');
for a = 1:2*t-1
    fprintf(fid, 'syndrome%u, ', a);
end
fprintf(fid, 'codeword, message);\n');
fprintf(fid, '  input [%u:0] received;\n', n-1);
for a = 1:2*t-1
    fprintf(fid, '  input [%u:0] syndrome%u;\n', m-1, a);
end
fprintf(fid, '  output [%u:0] codeword;\n', n-1);
fprintf(fid, '  output [%u:0] message;\n\n', k-1);
for a = 0:t
    fprintf(fid, '  wire [%u:0] locator%u;\n', m-1, a);
end
fprintf(fid, '\n');
fprintf(fid, '  bch_ibm ibm(');
for a = 1:2*t-1
    fprintf(fid, 'syndrome%u, ', a);
end
for a = 0:t
    if a
        fprintf(fid, ', ');
    end
    fprintf(fid, 'locator%u', a);
end
fprintf(fid, ');\n');
fprintf(fid, '  bch_chien chien(received');
for a = 0:t
    fprintf(fid, ', locator%u', a);
end
fprintf(fid, ', codeword);\n\n');
fprintf(fid, '  assign message = codeword[%u:%u];\n', n-1, n-k);
fprintf(fid, 'endmodule\n');
fclose(fid);

% Write the Verilog for the decoder
fid = fopen('bch_decoder.v', 'w');
fprintf(fid, 'module bch_decoder(received, codeword, message, ');
fprintf(fid, 'detection);\n');
fprintf(fid, '  input [%u:0] received;\n', n-1);
fprintf(fid, '  output [%u:0] codeword;\n', n-1);
fprintf(fid, '  output [%u:0] message;\n', k-1);
fprintf(fid, '  output detection;\n\n');
for a = 1:2*t-1
    fprintf(fid, '  wire [%u:0] syndrome%u;\n', m-1, a);
end
fprintf(fid, '\n');
fprintf(fid, '  bch_detector detector(received, ');
for a = 1:2*t-1
    fprintf(fid, 'syndrome%u, ', a);
end
fprintf(fid, 'detection);\n');
fprintf(fid, '  bch_corrector corrector(received, ');
for a = 1:2*t-1
    fprintf(fid, 'syndrome%u, ', a);
end
fprintf(fid, 'codeword, message);\n\n');
fprintf(fid, '  assign message = codeword[%u:%u];\n', n-1, n-k);
fprintf(fid, 'endmodule\n');
fclose(fid);

if make_testbench
    % Write the Verilog for the encoder testbench
    rs = RandStream('mt19937ar');
    fid = fopen('bch_encoder_tb.v', 'w');
    fprintf(fid, 'module encoder_tb;\n');
    fprintf(fid, '  reg [%u:0] message;\n', k-1);
    fprintf(fid, '  wire [%u:0] codeword;\n', n-1);
    fprintf(fid, '  integer failures;\n\n');
    fprintf(fid, '  bch_encoder encoder(message, codeword);\n\n');
    fprintf(fid, '  initial\n');
    fprintf(fid, '  begin\n');
    fprintf(fid, '    failures = 0;\n\n');
    for a = 1:20
        m = gf([zeros(1,l) randi(rs, [0 1], 1, k)]);
        u = bchenc(m, n+l, k+l);
        m = m.x;
        m = num2str(m(:,1+l:k+l), '%u');
        u = u.x;
        u = num2str(u(:,1+l:n+l), '%u');
        fprintf(fid, '    message <= %u''b%s;\n', k, m);
        fprintf(fid, '    #1;\n');
        fprintf(fid, '    if (codeword != %u''b%s)\n', n, u);
        fprintf(fid, '      failures = failures+1;\n\n');
    end
    fprintf(fid, '    $display("\\n==============================");\n');
    fprintf(fid, '    if (!failures)\n');
    fprintf(fid, '      $display("\\nAll encoding tests passed\\n");\n');
    fprintf(fid, '    else\n');
    fprintf(fid, ...
        '      $display("Failed %%d encoding test(s)", failures);\n');
    fprintf(fid, '    $display("==============================\\n");\n');
    fprintf(fid, '  end\n');
    fprintf(fid, 'endmodule\n');
    fclose(fid);
    
    % Write the Verilog for the decoder testbench
    fid = fopen('bch_decoder_tb.v', 'w');
    fprintf(fid, 'module decoder_tb;\n');
    fprintf(fid, '  reg [%u:0] received;\n', n-1);
    fprintf(fid, '  wire [%u:0] codeword;\n', n-1);
    fprintf(fid, '  wire [%u:0] message;\n', k-1);
    fprintf(fid, '  wire detection;\n');
    fprintf(fid, '  integer failures;\n\n');
    fprintf(fid, '  bch_decoder decoder(received, codeword, message, ');
    fprintf(fid, 'detection);\n\n');
    fprintf(fid, '  initial\n');
    fprintf(fid, '  begin\n');
    fprintf(fid, '    failures = 0;\n\n');
    for a = 0:t+1
        for b = 1:20
            m = gf([zeros(1,l) randi(rs, [0 1], 1, k)]);
            u = bchenc(m, n+l, k+l);
            m = m.x;
            m = num2str(m(:,1+l:k+l), '%u');
            u = u(:,1+l:n+l);
            e = gf(zeros(1,n));
            for c = 1:a
                d = randi(rs, [1 n]);
                while e(d) == gf(1)
                    d = randi(rs, [1 n]);
                end
                e(d) = gf(1);
            end
            x = u+e;
            x = num2str(x.x, '%u');
            fprintf(fid, '    received <= %u''b%s;\n', n, x);
            fprintf(fid, '    #1;\n');
            if ~a
                fprintf(fid, '    if (message != %u''b%s || ', k, m);
                fprintf(fid, 'detection)\n');
            elseif a < t+1
                fprintf(fid, '    if (message != %u''b%s || ', k, m);
                fprintf(fid, '!detection)\n');
            else
                fprintf(fid, '    if (!detection)\n');
            end
            fprintf(fid, '      failures = failures+1;\n\n');
        end
    end
    fprintf(fid, '    $display("\\n==============================");\n');
    fprintf(fid, '    if (!failures)\n');
    fprintf(fid, '      $display("\\nAll decoding tests passed\\n");\n');
    fprintf(fid, '    else\n');
    fprintf(fid, ...
        '      $display("Failed %%d decoding test(s)", failures);\n');
    fprintf(fid, '    $display("==============================\\n");\n');
    fprintf(fid, '  end\n');
    fprintf(fid, 'endmodule\n');
    fclose(fid);
end
end