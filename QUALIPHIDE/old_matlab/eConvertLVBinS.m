function [b] = eConvertLVBinS(fname)
%
%
    fid = fopen(fname,'r');
    a = fread(fid, 2, 'double', 0, 'b');
    Ncoeffs = a(1);
    Nchs = a(2);
    a = fread(fid, [ Ncoeffs, Nchs ], 'double', 0, 'b');
    dt = 1./a(1,:);
    offsets = a(2,:);
    gains = a(3:end,:);    
    a = fread(fid, [ Nchs, inf ], 'int16', 0, 'b');
    b = repmat(offsets', [1, size(a,2)]);
    for nn = 1:size(gains,1)
        b = b + repmat(gains(nn,:)',[1, size(a,2)]).*a.^nn;                    
    end
    b = [(0:size(b,2)-1)'*dt(1,1),b'];
		fclose(fid);
end