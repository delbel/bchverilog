function [ ] = hammingverilog(k, make_testbench)
%hammingverilog Creates Verilog codec for a shortened systematic Hamming code

% Find the code parameters
m = 1;
while 2^m-m-1 < k
    m = m+1;
end
l = 2^m-m-1-k;
n = 2^m-1-l;

% Write the Verilog for the encoder
[H, G] = hammgen(m);
G = G(1:k, 1:n);
fid = fopen('hamming_encoder.v', 'w');
fprintf(fid, 'module hamming_encoder(message, codeword);\n');
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

% Write the Verilog for the detector
H = H(:, 1:n);
fid = fopen('hamming_detector.v', 'w');
fprintf(fid, 'module hamming_detector(received, syndrome, detection);\n');
fprintf(fid, '  input [%u:0] received;\n', n-1);
fprintf(fid, '  output [%u:0] syndrome;\n', n-k-1);
fprintf(fid, '  output detection;\n\n');
for a = 1:n-k
    fprintf(fid, '  assign syndrome[%u] = ', a-1);
    if ~sum(H(a, :))
        fprintf(fid, '1''b0\n');
    else
        c = 0;
        for b = 1:n
            if H(a, b)
                if c
                    fprintf(fid, ' ^ ');
                end
                c = c+1;
                fprintf(fid, 'received[%u]', b-1);
            end
        end
        fprintf(fid, ';\n');
    end
end
fprintf(fid, '  assign detection = (| syndrome);\n');
fprintf(fid, 'endmodule\n');
fclose(fid);

% Write the Verilog for the corrector
fid = fopen('hamming_corrector.v', 'w');
fprintf(fid, 'module hamming_corrector(received, syndrome, codeword, ');
fprintf(fid, 'message);\n');
fprintf(fid, '  input [%u:0] received;\n', n-1);
fprintf(fid, '  input [%u:0] syndrome;\n', n-k-1);
fprintf(fid, '  output [%u:0] codeword;\n', n-1);
fprintf(fid, '  output [%u:0] message;\n\n', k-1);
for a = 1:n
    fprintf(fid, '  assign codeword[%u] = ', a-1);
    fprintf(fid, 'received[%u] ^ (~| (syndrome ^ ', a-1);
    fprintf(fid, '%u''b%s));\n', n-k, num2str(flipud(H(:, a)), '%u'));
end
fprintf(fid, '  assign message = codeword[%u:%u];\n', n-1, n-k);
fprintf(fid, 'endmodule\n');
fclose(fid);

% Write the Verilog for the decoder
fid = fopen('hamming_decoder.v', 'w');
fprintf(fid, 'module hamming_decoder(received, codeword, message, ');
fprintf(fid, 'detection);\n');
fprintf(fid, '  input [%u:0] received;\n', n-1);
fprintf(fid, '  output [%u:0] codeword;\n', n-1);
fprintf(fid, '  output [%u:0] message;\n', k-1);
fprintf(fid, '  output detection;\n\n');
fprintf(fid, '  wire [%u:0] syndrome;\n\n', n-k-1);
fprintf(fid, 'hamming_detector detector(received, syndrome, ');
fprintf(fid, 'detection);\n');
fprintf(fid, 'hamming_corrector corrector(received, syndrome, codeword, ');
fprintf(fid, 'message);\n');
fprintf(fid, 'endmodule\n');
fclose(fid);

if make_testbench
    % Write the Verilog for the encoder testbench
    rs = RandStream('mt19937ar');
    fid = fopen('hamming_encoder_tb.v', 'w');
    fprintf(fid, 'module encoder_tb;\n');
    fprintf(fid, '  reg [%u:0] message;\n', k-1);
    fprintf(fid, '  wire [%u:0] codeword;\n', n-1);
    fprintf(fid, '  integer failures;\n\n');
    fprintf(fid, '  hamming_encoder encoder(message, codeword);\n\n');
    fprintf(fid, '  initial\n');
    fprintf(fid, '  begin\n');
    fprintf(fid, '    failures = 0;\n\n');
    for a = 1:20
        m = randi(rs, [0 1], 1, k);
        u = encode(m, n, k, 'linear/binary', flipud(fliplr(G)));
        m = num2str(m, '%u');
        u = num2str(u, '%u');
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
    fid = fopen('hamming_decoder_tb.v', 'w');
    fprintf(fid, 'module decoder_tb;\n');
    fprintf(fid, '  reg [%u:0] received;\n', n-1);
    fprintf(fid, '  wire [%u:0] codeword;\n', n-1);
    fprintf(fid, '  wire [%u:0] message;\n', k-1);
    fprintf(fid, '  wire detection;\n');
    fprintf(fid, '  integer failures;\n\n');
    fprintf(fid, '  hamming_decoder decoder(received, codeword, ');
    fprintf(fid, 'message, detection);\n\n');
    fprintf(fid, '  initial\n');
    fprintf(fid, '  begin\n');
    fprintf(fid, '    failures = 0;\n\n');
    for a = 0:2
        for b = 1:20
            m = randi(rs, [0 1], 1, k);
            u = gf(encode(m, n, k, 'linear/binary', flipud(fliplr(G))));
            m = num2str(m, '%u');
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
            elseif a < 2
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
