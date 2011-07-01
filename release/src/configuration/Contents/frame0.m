function Content = frame0(info)

name = 'frame0';

Content = content_template(name);
Content.wordlen = {[3],[4],[5],[6]};
Content.stepping = 3;

Content.offset = mod(3-mod(Content.wordlen{1},3),3);
 
