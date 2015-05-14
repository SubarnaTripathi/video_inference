%%convert_image

prefix = '001TP_';

for i = 63:100
    a = sprintf('%03d', i);
    ifilename = sprintf('%s%s.png',prefix,a);
    ofilename = sprintf('%s%s.jpg',prefix,a);
    
    im = imread(ifilename,'png');
    imwrite(im, ofilename, 'jpg');
end