def getCord(a,b,c):
	return float((c-a))/float((b-a))

def getAsympts(f11,f12,f21,f22):
	x = (f11-f12)/(f22+f11-f12-f21)
	y = (f11-f21)/(f22+f11-f12-f21)
	return x,y

if __name__ == "__main__":
	print(getCord(1,-2,0))
	print(getAsympts(-2.,4.,1.,-2.))