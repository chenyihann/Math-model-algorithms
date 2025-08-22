function Y=jieduan(X,l,w,a,b,L)
%  * brief   计算截断效率
%  * input   X为定日镜中心坐标(可为矩阵),l,w为定日镜长宽(可为矩阵),a,b为太阳
%  * 高度角和方位角,L为集热器中心坐标
%  * output  Y为截断效率(可为向量)

N=size(X,1);%计算总数
Y=zeros(N,1);%初始化Y

x=L(1);
y=L(2);
z=L(3);%L坐标

Sr=-[cos(a)*sin(b);cos(a)*cos(b);sin(a)];
%中心入射光线单位向量

T1=[cos(b),-sin(b)*sin(a),sin(b)*cos(a);
    sin(b),cos(b)*sin(a),-cos(b)*cos(a);
    0,cos(a),sin(a)];%光线坐标系到地面坐标系转换矩阵

for i=1:N
    Sf=[x-X(i,1);y-X(i,2);z-X(i,3)]./sqrt((X(i,1)-x)^2+(X(i,2)-y)^2+(X(i,3)-z)^2);
    %中心反射光线单位向量

    H=(Sf-Sr)./norm(Sf-Sr);%法向量
    A=acos(H(3));%镜面竖直倾斜角
    B=atan(H(1)/H(2));
    if(H(1)>=0 && H(2)>=0)
        B=B;
    elseif(H(1)>=0 && H(2)<0)
        B=B+pi;
    elseif(H(1)<0 && H(2)<0)
        B=B+pi;
    else
        B=B+2*pi;
    end%镜面方位角(即从y轴顺时针旋转到对应角的最小角度)

    T2=[cos(B),-sin(B)*sin(A),sin(B)*cos(A);
    sin(B),cos(B)*sin(A),-cos(B)*cos(A);
    0,cos(A),sin(A)];%光线坐标系到地面坐标系转换矩阵

    num=0;%记数
   for k=1:10000 %蒙特卡洛随机10000次
       p=acos(1-rand(1)*(1-cos(4.65*1e-3)));%随机极角
       q=2*pi*rand(1);%随机方位角

       M=-T1*[sin(p)*cos(q);sin(p)*sin(q);cos(p)];%随机入射光线单位向量
       N=-2*dot(M,H)*H+M;%对应反射光线单位向量

       C=T2*[-0.5*l(i)+rand(1)*l(i);-0.5*w(i)+rand(1)*w(i);0];%随机落点

       ap=N(1);bp=N(2);cp=N(3);
       xp=C(1);yp=C(2);zp=C(3);

       delta=(ap*xp+bp*yp)^2-(ap^2+bp^2)*(xp^2+yp^2-49);
       if(delta>0)
            answer=zp-cp*((ap*xp+bp*yp)+sqrt(delta))/(ap^2+bp^2);%交点z坐标
            if(z-4<=answer && answer<=z+4)%吸收判定条件
                num=num+1;
            end
       end
   end

   Y(i)=num/10000;
end
end