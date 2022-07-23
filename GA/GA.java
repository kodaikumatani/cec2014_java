import java.lang.Math;
import java.util.Random;
import java.util.Map;
import java.util.Set;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
/* csv */
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.IOException;

public class GA {
	public static int f = 5;	// 0, 1, 2, 3, ..., 22
	public static int dim = 10;	// 10 (or 30, 50, 100)
	public static double MIN = -100.0;
	public static double MAX = 100.0;
	
	public static void main(String[] args) {
		//System.out.println(calulation(100000, 50, 15, 12));
		
		experiment();
		
		/*
		double dist = Landscape.distance(f, x);
		System.out.println("Distance = " + dist);
		double[] fval = new double[Landscape.numpoints+1];
		fval = Landscape.evaluate(f, vector);
		for(int i=0; i<=Landscape.numpoints; i++){
			System.out.println(fval[i]);
		}
		*/
	}
	
	private static void experiment()	{
		int trial = 51;
		int eval_times = 100000;

		List<Integer> ramda = Arrays.asList(20,50,80,125,160);
		List<Double> mu = Arrays.asList(0.2,0.4,0.6,0.8);
		List<Double> nu = Arrays.asList(0.2,0.4,0.6,0.8);
		Map<String, String[]> type = new TreeMap<>();
		type.put("uniform-1_0", new String[]{"uniform","1.0"});
		type.put("uniform-2_0", new String[]{"uniform","2.0"});
		type.put("uniform-3_0", new String[]{"uniform","3.0"});
		type.put("gaussian-0_4", new String[]{"gaussian", "0.4"});
		type.put("gaussian-0_8", new String[]{"gaussian", "0.8"});
		type.put("gaussian-1_2", new String[]{"gaussian", "1.2"});
		
		type.forEach((key, value) -> {
			System.out.println(key);
			int count = 1;
			List<Double[]> data = new ArrayList<>();
			for (int i : ramda) {
				for (Double j : mu) {
					for (Double k : nu) {
						int discovered = 0;
						Double[] eval = new Double[trial];
						for (int l = 0; l < trial; l++) {
							eval[l] = calulation(eval_times, i, (int)(i*j), (int)(i*k), value);
							if (eval[l] < 520) discovered++;
						}
						Double[] cpy = Arrays.copyOf(eval,trial);
						Arrays.sort(cpy);
						double median = cpy[cpy.length / 2];
						
						Double[] line = new Double[trial+4];
						System.arraycopy(new Double[]{(double)i, j, k, median}, 0, line, 0, 4);
						System.arraycopy(eval, 0, line, 4, trial);
						//System.arraycopy(コピー元配列, コピー元配列のコピー開始位置, コピー先配列, コピー先配列の開始位置, コピーの個数)
						data.add(line);
						progressBar(count++, ramda.size() * mu.size() * nu.size());
					}
				}
			}
			exportCsv(data, key);
		});
	}
	

	public static double calulation(int eval_times, int ramda, int mu, int nu, String[] status)	{
		int max_gene = (int) eval_times / ramda;
		Map<Double, double[]> children = getFirstChildren(ramda);
		Map<Double, double[]> elites = selection(nu, children);
		Map<Double, double[]> parents = selection(mu, children);
		for (int i = 1; i < max_gene; i++)	{
			children = getChildren(ramda, parents, elites, status);
			elites = selection(nu, concatenate(elites, children));
			parents = selection(mu, children);
		}
		return elites.keySet().toArray(new Double[elites.size()])[0];
	}

	private static Map<Double, double[]> getFirstChildren(int ramda) {
		Map<Double, double[]> children = new TreeMap<>();
		for (int i = 0; i < ramda; i++) {
			double[] vector = new Random().doubles(dim, MIN, MAX).toArray();
			double value = Simulator.evaluate(f, vector);
			children.put(value, vector);
		}
		return children;
	}

	private static Map<Double, double[]> selection(int num, Map<Double, double[]> map) {
		Map<Double, double[]> selected = new TreeMap<>();
		Set<Double> keys = map.keySet();
		for (int i = 0; i < num; i++) {
			Double key = keys.toArray(new Double[0])[i];
			//System.out.println(key + " => " + map.get(key));
			selected.put(key,map.get(key));
		}
		return selected;
	}

	private static Map<Double, double[]> getChildren(int ramda, Map<Double, double[]> parents, Map<Double, double[]> elites, String[] status) {		
		Map<Double, double[]> children = new HashMap<>();
		Map<Double, double[]> concat = concatenate(parents, elites);
		for (int i = 0; i < ramda; i++)	{
			//double[] p_vector = BrendCrossover(concat);
			double[] p_vector = TriangleCrossover(concat, status);
			//カットオフを行う
			double[] vector = new double[dim];
 			for(int j = 0; j < dim; j++)	{
				double coordinate = p_vector[j];
				coordinate = coordinate > MAX ? MAX : coordinate;
				coordinate = coordinate < MIN ? MIN : coordinate;
				vector[j] = coordinate;
			}
			vector = mutation(vector);
			double value = Simulator.evaluate(f, vector);
			children.put(value, vector);
		}
		//System.out.println(randomValue.getClass().getSimpleName());
		return children;
	}

	private static double[] BrendCrossover(Map<Double, double[]> map)	{
		double alpha = 0.5;
		// ２点を選択する
		List<double[]> vertex = new ArrayList<>();
		Random generator = new Random();
		List<Double> keys = new ArrayList<>(map.keySet());
		for (int i = 0; i < 2; i++)	{
			int random = generator.nextInt(keys.size());
			vertex.add(map.get(keys.get(random)));
			keys.remove(random);
		}
		
		double[] list = new double[dim];
		for (int i = 0; i < dim; i++)	{
			double A = vertex.get(0)[i];
			double B = vertex.get(1)[i];
			double midpoint = (A + B) / 2;
			double L = A > B ? Math.abs(A-B) : Math.abs(B-A);
			A = midpoint + (L/2+L*alpha);
			B = midpoint - (L/2+L*alpha);
			list[i] = B + L * generator.nextDouble();
			//System.out.println(L/2+L*alpha);
		}
		return list;
	}

	private static double[] TriangleCrossover(Map<Double, double[]> map, String[] status)	{
		// ３点を選択する
		List<double[]> vertex = new ArrayList<>();
		List<Double> keys = new ArrayList<>(map.keySet());
		for (int i = 0; i < 3; i++)	{
			int random = new Random().nextInt(keys.size());
			vertex.add(map.get(keys.get(random)));
			keys.remove(random);
		}

		double value = Double.parseDouble(status[1]);
		
		double k = 0;
		double l = 0;
		double m = 0;
		if (status[0] == "uniform")	{
			k = new Random().nextDouble(value);
			l = new Random().nextDouble(value);
			m = new Random().nextDouble(value);
		} else if (status[0] == "gaussian") {
			k = new Random().nextGaussian() * value;
			l = new Random().nextGaussian() * value;
			m = new Random().nextGaussian() * value;
		} else {
			System.out.println("Error:Not select random type");
			System.exit(0);
		}

		double[] list = new double[dim];
		for (int i = 0; i < dim; i++) {
			double OG = (  vertex.get(0)[i] + vertex.get(1)[i] + vertex.get(2)[i]) / 3;
			double AG = (2*vertex.get(0)[i] - vertex.get(1)[i] - vertex.get(2)[i]) / 3;
			double BG = (2*vertex.get(1)[i] - vertex.get(0)[i] - vertex.get(2)[i]) / 3;
			double CG = (2*vertex.get(2)[i] - vertex.get(0)[i] - vertex.get(1)[i]) / 3;
			list[i] = k*AG + l*BG + m*CG + OG;
		}
		return list;
	}

	private static double[] mutation(double[] vector)	{
		Random generator = new Random();
		double[] list = new double[dim];
 		for (int i = 0; i < vector.length; i++)	{
			int p = generator.nextInt(dim * 10);
			if (p % dim == 0)	{
				list[i] = MIN + (MAX - MIN) * generator.nextDouble();
			} else {
				list[i] = vector[i];
			}
		}
		return list;
	}

	private static Map<Double, double[]> concatenate(Map<Double, double[]> map_A, Map<Double, double[]> map_B) {
		Map<Double, double[]> map = new TreeMap<>();
		map_A.forEach((key, value) -> {
			map.put(key, value);
		});
		map_B.forEach((key, value) -> {
			map.put(key, value);
		});
		return map;
	}

	private static void exportCsv(List<Double[]> data, String filename)	{
        try {
            // 出力ファイルの作成
            FileWriter fw = new FileWriter("data/imp1/func_"+Integer.valueOf(f).toString()+"/"+filename+".csv", false);
            // PrintWriterクラスのオブジェクトを生成
            PrintWriter pw = new PrintWriter(new BufferedWriter(fw));
 
            // データを書き込む
            for(Double[] line : data) {
				for (int i = 0; i < line.length; i++)	{
					pw.print(line[i]);
					pw.print(",");
				}
				pw.println();
            }
            // ファイルを閉じる
            pw.close();
 
            // 出力確認用のメッセージ
            System.out.println("csvファイルを出力しました");
 
        // 出力に失敗したときの処理
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

	public static void progressBar(int i, int trial) {
        double d = (double)i/trial;
		String p = String.format("%6.2f", d * 100);
		double p_cnt = (d * 100) / 10;
		p_cnt *= 2; // 100%で20個表示
		StringBuffer buf = new StringBuffer();
		buf.append(p);
		buf.append("% [");
		for (int j = 0; j < 20; j++) {
			if (j < (int)p_cnt) {
				buf.append("=");
			} else {
				buf.append(" ");
			}
		}
		buf.append("]");
		buf.append(Integer.toString(i));
		if (i == trial) {
			buf.append("\n");
		} else {
			buf.append("\r");
		}
		System.out.print(buf.toString());
	}
}