
rapcf_density_fun <- function( data, v, m = NULL, student = FALSE, transform = exp ) {
  # assume variances are normal 
  
  sigma = transform(v)
  
  if( !student ) {
    return( mnormt::dmnorm( data, varcov = sigma ))
  }
  # otherwise use T dist
  else {
    # mnormt::dmt for multivariate T dist
    # TODO: translate
    # [d d N] =size(Sigma);    
    # modelDensity=zeros(N,1);
    # ie if there is a parameter 
    # nu=exp(m(:,end))+2; #%log(nu-2)~N(.,.)
    # for i=1:N
    # covAdj= (nu(i)-2)./nu(i);
    # modelDensity(i,:)=mvtpdfln(data,Sigma(:,:,i).*covAdj,nu(i),[]);
  }
} 


