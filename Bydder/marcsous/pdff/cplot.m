function cplot(x,y,z)

h = ishold;

if nargin==1
    xx = 1:length(x);
    yy = x;
    ropts = '';
    iopts = '--';
end
if nargin==3
    xx = x;
    yy = y;
    ropts = z;
    iopts = z;
end
if nargin==2
    if ischar(y)
        xx = 1:length(x);
        yy = x;
        ropts = y;
        iopts = y;
    else
        xx = x;
        yy = y;
        ropts = '';
        iopts = '--r';
    end
end

% force red/blue color scheme
blue = [0.0000 0.4470 0.7410];
red  = [0.8500 0.3250 0.0980];

plot(xx,real(yy),ropts,'Color',blue);
hold on
plot(xx,imag(yy),iopts,'Color',red);

if h==0
    hold off
end