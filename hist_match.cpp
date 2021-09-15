Mat histogram_match(Mat image, Mat image3){
	string str2="/home/paras/ME3.jpg";
	string str="/home/paras/ME5.jpg";
	Mat image2;
	double intensities[256]={0};
	double pi;
	double intensities2[256]={0};
	double pi2;
	double intensities3[256]={-1};
	int max1=0, min1=255, max2=0, min2=255, tmp;
	long long sum=0;
	image2 = Mat::zeros(image.rows, image.cols, CV_8UC1);
	int size=image.rows*image.cols;
	int size2=image3.rows*image3.cols;
	for(int r=0; r<image.rows; r++){
		for(int c=0; c<image.cols; c++){
			intensities[image.at<uint8_t>(r,c)]+=1;
			if(max2<image.at<uint8_t>(r,c)){
				max2=image.at<uint8_t>(r,c);
			}
			if(min2>image.at<uint8_t>(r,c)){
				min2=image.at<uint8_t>(r,c);
			}
		}
	}
	for(int r=0; r<image3.rows; r++){
		for(int c=0; c<image3.cols; c++){
			intensities2[image3.at<uint8_t>(r,c)]+=1;
		}
	}
	for(int i=0; i<256; i++){
		intensities[i]+=intensities[i-1];
		intensities2[i]+=intensities2[i-1];
	}
	for(int i=0; i<256; i++){
		pi=intensities[i];
		pi=(pi*255)/size;
		intensities[i]=pi;
		pi=(int)intensities[i];
		if(intensities[i]-pi>.5f){
			intensities[i]=pi+1;
		}
		else{
			intensities[i]=pi;
		}
		pi2=intensities2[i];
		pi2=(pi2*255)/size2;
		intensities2[i]=pi2;
		pi2=(int)intensities2[i];
		if(intensities2[i]-pi2>.5f){
			intensities2[i]=pi2+1;
		}
		else{
			intensities2[i]=pi2;
		}
		if(intensities[i]>=255){
			intensities[i]=255;
		}
		if(intensities2[i]>=255){
			intensities2[i]=255;
		}
	}
	for(int i=0; i<256; i++){
        int j = 0;
        do {
              intensities3[i] = j;
              j++;
        }while(intensities[i] > intensities2[j]);

	}

	for(int i=0; i<256; i++){
		if(max1<intensities3[i]){
			max1=intensities3[i];
		}
		if(min1>intensities3[i]&&intensities3[i]!=-1){
			min1=intensities3[i];
		}
	}
//cout<<"     "<<max<<" "<<min<<" "<<max2<<" "<<min2<<endl;
//for(int i=0; i<256; i++)
//	cout<<i<<" "<<intensities[i]<<" "<<intensities2[i]<<" "<<intensities3[i]<<endl;

	for(int r=0; r<image2.rows; r++){
		for(int c=0; c<image2.cols; c++){
			if(intensities3[image.at<uint8_t>(r,c)]==-1){
//				tmp=image.at<uint8_t>(r,c);
				tmp=0;
			}
			else{
				tmp=intensities3[(int)image.at<uint8_t>(r,c)];
//				if(image.at<uint8_t>(r,c)==255){
//					cout<<"a="<<(uint8_t)image.at<uint8_t>(r,c)<<"c="<<tmp<<endl;
//				}
			}
//			tmp=((tmp-min1)*(max2-min2)/(max1-min1));
			image2.at<uint8_t>(r,c)=tmp;
			sum+=min((int)tmp,(int)intensities2[(int)image.at<uint8_t>(r,c)]);
//			sum+=tmp;
		}
	}
//	cout<<"sum="<<sum<<endl;
	return image2;
}
