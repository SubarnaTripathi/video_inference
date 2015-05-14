%%%%
function read_msh()

[width, height, band, img] = read_show('001TP_085.msh');
disp('read done');
write_im(width, height, band, img,'001TP_085_dup.msh');
disp('write done');
[width, height, band, img1] = read_show('001TP_085_dup.msh');
disp('second read');
end

function [width, height, band, img] = read_show(name)
fid = fopen(name, 'r');
width = fread(fid, 1,'uint32'); %'integer*4';
height = fread(fid, 1,'uint32'); %'integer*4';
band = fread(fid, 1,'uint32'); %'integer*4';

A = fread(fid, width*height, 'uint32'); %'integer*4';
fclose(fid);

B = reshape(A, width, height);
img = mat2gray(B');
imgUpsideDown = flipdim(img,1);
figure(1), imshow(imgUpsideDown);
end

function write_im(width, height, band, img,write_name)

fid = fopen(write_name, 'w');
fwrite(fid, width,'uint32'); %'integer*4';
fwrite(fid, height,'uint32'); %'integer*4';
fwrite(fid, band,'uint32'); %'integer*4';

img1 = uint32(255*img);

%img2 = img1;
figure(2), imshow(uint8(img1));

fwrite(fid, img1', 'uint32'); %'integer*4';
fclose(fid);
end
