function llh=xyz2llh(xyz,datum)
%XYZ2LLH   Calculates longitude, latitude, and height from global cartesisan coordinates.
%   LLH=xyz2llh(XYZ,DATUM) calculates longitude (deg), latitude (deg), and
%   height (m) on the ellipsoid specified by DATUM from the global cartestian
%   coordinates given in the 3xn (n = number of coordinate triples) matrix XYZ.
%   DATUM can either be a vector the first two elements of which give da and df,
%   or it can be a string containing the name of a datum that is resolved
%   by the function DATUMS function.
%
%   Note that longitude occupies the first row of LLH.
%
%   See DATUMS for more information on datum parameters.

%-------------------------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Aug 20, 2001  Peter Cervelli        Standardized code
%   Unknown       Peter Cervelli		Original Code
%
%-------------------------------------------------------------------------------

%Check input arguments

    if nargin==1
        da=0;
        df=0;
    else
        if strcmp(class(datum),'char')
            datum=datums(datum);
            if any(isnan(datum))
                error('Could not resolve datum name.')
            end
        end
        da=datum(1);
        df=datum(2);
    end
	
	if size(xyz,1)~=3
        error('Input xyz MUST be 3xn.')
	end

%Set constants

    a=6378137 - da;
    f=1/298.2572235630 - df;		
    b=(1-f)*a;
	e2=2*f - f^2;
	E2=(a^2-b^2)/(b^2);
   
%Calculate longitude, latitude, and height

	p=sqrt(xyz(1,:).^2 + xyz(2,:).^2);
    llh(1,:)=atan2(xyz(2,:),xyz(1,:));
    theta=atan((xyz(3,:)*a)./(p*b));
    llh(2,:)=atan((xyz(3,:)+E2*b*sin(theta).^3)./(p-e2*a*cos(theta).^3) );
    N=a./sqrt(1-e2*sin(llh(2,:)).^2);
    llh(3,:)=p./cos(llh(2,:))-N;

%Convert to degrees

    llh(1:2,:)=llh(1:2,:)*57.295779513082323;