function out = fetch_para(n,opt)

if opt==2
    syms m1;syms m2;
    syms m1rx1;syms m1ry1;syms m1rz1;
    syms m2rx2;syms m2ry2;syms m2rz2;

    
    syms xx1;syms xy1;syms xz1;syms yy1;syms yz1;syms zz1;
    syms xx2;syms xy2;syms xz2;syms yy2;syms yz2;syms zz2;

    m=[m1,m2];
    r=[m1rx1, 0,  0;
      m2rx2,  0,  0;];
    I=[zz1/2,0,0;
       0,zz1/2,0;
       0,0,-zz1/2;

       zz2/2,0,0;
       0,zz2/2,0;
       0,0,-zz2/2;
       ];

    % 
    mi=m(n);
    ri=r(n,:).';
    Ii=I(3*n-2:3*n,:);
    out=[Ii,ri;ri.',mi];
    
else
    syms m1;syms m2;syms m3;syms m4;syms m5;syms m6;syms m7;
    syms m1rx1;syms m1ry1;syms m1rz1;
    syms m2rx2;syms m2ry2;syms m2rz2;
    syms m3rx3;syms m3ry3;syms m3rz3;
    syms m4rx4;syms m4ry4;syms m4rz4;
    syms m5rx5;syms m5ry5;syms m5rz5;
    syms m6rx6;syms m6ry6;syms m6rz6;
    syms m7rx7;syms m7ry7;syms m7rz7;
    %
    syms xx1;syms xy1;syms xz1;syms yy1;syms yz1;syms zz1;
    syms xx2;syms xy2;syms xz2;syms yy2;syms yz2;syms zz2;
    syms xx3;syms xy3;syms xz3;syms yy3;syms yz3;syms zz3;
    syms xx4;syms xy4;syms xz4;syms yy4;syms yz4;syms zz4;
    syms xx5;syms xy5;syms xz5;syms yy5;syms yz5;syms zz5;
    syms xx6;syms xy6;syms xz6;syms yy6;syms yz6;syms zz6;
    syms xx7;syms xy7;syms xz7;syms yy7;syms yz7;syms zz7;
    m=[m1,m2,m3,m4,m5,m6,m7];
    r=[m1rx1, m1ry1,  m1rz1;
      m2rx2,  m2ry2,  m2rz2;
      m3rx3,  m3ry3,  m3rz3;
      m4rx4,  m4ry4,  m4rz4;
      m5rx5,  m5ry5,  m5rz5;
      m6rx6,  m6ry6,  m6rz6;
      m7rx7,  m7ry7,  m7rz7;];
    I=[(-xx1+yy1+zz1)/2,xy1,xz1;
       xy1,(xx1-yy1+zz1)/2,yz1;
       xz1,yz1,(xx1+yy1-zz1)/2;

       (-xx2+yy2+zz2)/2,xy2,xz2;
       xy2,(xx2-yy2+zz2)/2,yz2;
       xz2,yz2,(xx2+yy2-zz2)/2;

       (-xx3+yy3+zz3)/2,xy3,xz3;
       xy3,(xx3-yy3+zz3)/2,yz3;
       xz3,yz3,(xx3+yy3-zz3)/2;

       (-xx4+yy4+zz4)/2,xy4,xz4;
       xy4,(xx4-yy4+zz4)/2,yz4;
       xz4,yz4,(xx4+yy4-zz4)/2;

       (-xx5+yy5+zz5)/2,xy5,xz5;
       xy5,(xx5-yy5+zz5)/2,yz5;
       xz5,yz5,(xx5+yy5-zz5)/2;

      (-xx6+yy6+zz6)/2,xy6,xz6;
       xy6,(xx6-yy6+zz6)/2,yz6;
       xz6,yz6,(xx6+yy6-zz6)/2;

      (-xx7+yy7+zz7)/2,xy7,xz7;
       xy7,(xx7-yy7+zz7)/2,yz7;
       xz7,yz7,(xx7+yy7-zz7)/2;];

    
    mi=m(n);
    ri=r(n, :).';
    Ii=I(3*n-2:3*n, :);
    out=[Ii,ri; ri.', mi];
end
    
end

