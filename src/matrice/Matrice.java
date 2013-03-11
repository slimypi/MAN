/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package matrice;

/**
 *
 * @author Slimane
 */
public class Matrice {

   private final int M;             // number of rows
    private final int N;             // number of columns
    private final double[][] data;   // M-by-N array

    // create M-by-N matrix of 0's
    public Matrice(int M, int N) {
        this.M = M;
        this.N = N;
        data = new double[M][N];
        for(int i =0;i<M;i++)  //initialisation a zéro
        {
            for(int j=0;j<N;j++)
            {
                this.data[i][j]=0.0;
            }
        }
    }

    
    
    // create matrix based on 2d array
    public Matrice(double[][] data) {
        M = data.length;
        N = data[0].length;
        this.data = new double[M][N];
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                    this.data[i][j] = data[i][j];
            }
        }
    }

    // copy constructor
    private Matrice(Matrice A) { this(A.data); }

    // create and return a random M-by-N matrix with values between 0 and 1
    public static Matrice random(int M, int N) {
        Matrice A = new Matrice(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                A.data[i][j] = Math.random();
        return A;
    }

    public void remplire(int i,int j, double n)
    {
        this.data[i][j]=n;
    }
    
    // create and return the N-by-N identity matrix
    public static Matrice identity(int N) {
        Matrice I = new Matrice(N, N);
        for (int i = 0; i < N; i++)
            I.data[i][i] = 1;
        return I;
    }

   
    // swap rows i and j
    private void swap(int i, int j) {
        double[] temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }

    // create and return the transpose of the invoking matrix
    public Matrice transpose() {
        Matrice A = new Matrice(N, M);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                A.data[j][i] = this.data[i][j];
        return A;
    }

    // return C = A + B
    public Matrice plus(Matrice B) {
        Matrice A = this;
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        Matrice C = new Matrice(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] + B.data[i][j];
        return C;
    }


    // return C = A - B
    public Matrice minus(Matrice B) {
        Matrice A = this;
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        Matrice C = new Matrice(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] - B.data[i][j];
        return C;
    }

    // does A = B exactly?
    public boolean eq(Matrice B) {
        Matrice A = this;
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                if (A.data[i][j] != B.data[i][j]) return false;
        return true;
    }

    // return C = A * B
    public Matrice times(Matrice B) {
        Matrice A = this;
        if (A.N != B.M) throw new RuntimeException("Illegal matrix dimensions.");
        Matrice C = new Matrice(A.M, B.N);
        for (int i = 0; i < C.M; i++)
            for (int j = 0; j < C.N; j++)
                for (int k = 0; k < A.N; k++)
                    C.data[i][j] += (A.data[i][k] * B.data[k][j]);
        return C;
    }


    // return x = A^-1 b, assuming A is square and has full rank
    public Matrice solve(Matrice rhs) {
        if (M != N || rhs.M != N || rhs.N != 1)
            throw new RuntimeException("Illegal matrix dimensions.");

        // create copies of the data
        Matrice A = new Matrice(this);
        Matrice b = new Matrice(rhs);

        // Gaussian elimination with partial pivoting
        for (int i = 0; i < N; i++) {

            // find pivot row and swap
            int max = i;
            for (int j = i + 1; j < N; j++)
                if (Math.abs(A.data[j][i]) > Math.abs(A.data[max][i]))
                    max = j;
            A.swap(i, max);
            b.swap(i, max);

            // singular
            if (A.data[i][i] == 0.0) throw new RuntimeException("Matrix is singular.");

            // pivot within b
            for (int j = i + 1; j < N; j++)
                b.data[j][0] -= b.data[i][0] * A.data[j][i] / A.data[i][i];

            // pivot within A
            for (int j = i + 1; j < N; j++) {
                double m = A.data[j][i] / A.data[i][i];
                for (int k = i+1; k < N; k++) {
                    A.data[j][k] -= A.data[i][k] * m;
                }
                A.data[j][i] = 0.0;
            }
        }

        // back substitution
        Matrice x = new Matrice(N, 1);
        for (int j = N - 1; j >= 0; j--) {
            double t = 0.0;
            for (int k = j + 1; k < N; k++)
                t += A.data[j][k] * x.data[k][0];
            x.data[j][0] = (b.data[j][0] - t) / A.data[j][j];
        }
        System.out.println("Sortie Marice Vecteur X: ");
        return x;
   
    }

    // print matrix to standard output
    public void show() {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) 
                System.out.printf("%9.4f ", data[i][j]);
            System.out.println();
        							}
    					}
    
    /*public void reduction_Gauss (Matrice b)
    {
    	
    	for (int k=0;k<b.N-1;k++)
    	{
    		for (int i=k+1;i<=N-1;i++ )
    		{
    			double c = this.data[i-1][k-1]/this.data[k-1][k-1];
    			b.data[i-1][0]=b.data[i-1][0] - c*b.data[k][0];
    			this.data[i-1][k-1]=0;
    		
    			for(int j=k+1;j<=b.N-1;j++)
    			{
    				this.data[i-1][j-1]=this.data[i-1][j-1] - c*this.data[k-1][j-1];
    				this.
    			}
    		}
    	}
    	
    }
}
    
    public Matrice[] reductionLU()
    {
    	Matrice tableau[]= new Matrice[2];
    	tableau[0]= new Matrice(N,M);
    	tableau[1]= new Matrice(N,M);
    	for (int i = 0; i < N; i++)
            tableau[0].data[i][i] = 1;
    	int i;
    	for ( i =1;i<N;i++)
    	{
    		for(int j=i;j<=N;j++)
    		{
    		
    			double s=0;
    			for(int k=1;k<=i-1;k++)
    			{
    				s=s+tableau[0].data[i][k]*tableau[1].data[k][j];
    			}
    		
    			tableau[1].data[i][j]=this.data[i][j]-s;
    		}
    		
    		for(int j = i+1;j<=this.N;j++)
    		{
    			
    			double ss=0;
    			for(int k=1;k<=i-1;k++)
    			{
    				
    				ss=ss+tableau[0].data[j][k]*tableau[1].data[k][i];
    			}
    			
    			
				tableau[0].data[j][i]= 1/tableau[1].data[i][i]*(this.data[j][i]- ss);
    		}
    	}
    	for(int j = i+1;j<=this.N;j++)
		{
			
			double sss=0;
			for(int k=1;k<=this.N-1;k++)
			{
				
				sss=sss+tableau[0].data[this.N][k]*tableau[1].data[k][this.N];
			}
			
    			tableau[1].data[this.N][this.N]=this.data[N][N]-sss;
    			
    }
    	return tableau;
    }*/
    
    public Matrice getL()
    { 
    	Matrice X=new Matrice(N,M);
    	double [][] L = X.data;
    	for(int i=0;i<M;i++)
    	{
    		for(int j=0;j<N;j++)
    		{
    			if(i>j)
    			{
    				L[i][j]=data[i][j];
    			}
    			else if (i==j)
    			{
    				L[i][j]=1.0;
    			}
    			else
    			{
    				L[i][j]=0.0;
    			}
    		}
    	}
        System.out.println("Sortie Marice Triangulaire Inférieure : ");
    	return X;
    }
    
    public Matrice getU()
    {
    	Matrice X=new Matrice(N,N);
    	double[][] U = X.data;
    	for(int i=0;i<N;i++)
    	{
    		for(int j=0;j<N;j++)
    		{
    			if(i<=j)
    			{
    				U[i][j]=data[i][j];
    			}
    			else
    			{
    				U[i][j]=0.0;
    			}
    		}
    	}
        System.out.println("Sortie Marice Triangulaire Supérieure : ");
    	return X;
    }
    
    public Matrice cholesky()
    {
        int m = this.M;
        Matrice X = new Matrice(m,m); 
        for(int i =0;i<m;i++)
        {
            for(int k =0;k<=i;k++)
            {
                double sum=0;
                for(int j=0;j<k;j++)
                {
                    sum+=X.data[i][j]*X.data[k][j];
                }
            if(i==k) X.data[i][i]=Math.sqrt(this.data[i][i]- sum);
            else X.data[i][k]=1.0/X.data[k][k] * (X.data[i][k] - sum);
            }
        }
        System.out.println("Sortie Marice Cholesky : ");
        return X;
    }
    public static void main(String[] args) {
        Matrice A = Matrice.random(5,5);
		A.show();
		Matrice b = Matrice.random(5, 1);
        //b.show();
        System.out.println();

        Matrice x = A.solve(b);
        x.show();
        //System.out.println("Solution: Matrice Triangulaire Supérieure");
        Matrice U,L;
        U=A.getU();
        U.show();
        L=A.getL();
        //System.out.println("Solution: Matrice Triangulaire Inférieure");
        L.show();
        //System.out.println("Cholesky ");
        Matrice C;
        C=A.cholesky();
        C.show();
        Matrice D=new Matrice(3,3);
        D.remplire(0,0,12);
        D.remplire(1,1,11);
        D.remplire(2, 2, 15);
        Matrice btest = new Matrice(3,1);
        btest.remplire(0,0,1);
        btest.remplire(1, 0, 1);
        btest.remplire(2,0, 1);
        btest.show();
        D.show();
        //Matrice XX=new Matrice(1,3);
       Matrice XX=D.solve(btest);
        XX.show();
        
        
     }
}
