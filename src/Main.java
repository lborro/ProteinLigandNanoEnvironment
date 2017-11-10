import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import br.embrapa.cnptia.cbi.sdl.utils.Utils;
import br.embrapa.cnptia.gpbc.plc.data.MaxContactsTable;
import br.embrapa.cnptia.gpbc.plc.descriptors.protein.ProteinNanoEnvironment;


public class Main {

	public static void main(String[] args) throws Exception {
		//System.setProperty("http.proxyHost", "proxy.cnptia.embrapa.br");
		//System.setProperty("http.proxyPort", "3128");

		ExecutorService pool = Executors.newFixedThreadPool(3);
		Set<Future<ProteinNanoEnvironment>> results = new HashSet<Future<ProteinNanoEnvironment>>();

		MaxContactsTable maxCon = new MaxContactsTable("data/maxcon.txt");

		BufferedReader br = new BufferedReader(new FileReader(args[0]));
		String line;
		while ((line = br.readLine()) != null) {
			String[] tk = line.split("\t");
			File protein = null, ligand = null;
			try{
				protein = new File(tk[0].trim());
				ligand = new File(tk[1].trim());
				Callable<ProteinNanoEnvironment> callable = new NanoEnvironmentCalculation(protein, ligand, maxCon, 5) ;
				Future<ProteinNanoEnvironment> future = pool.submit(callable);
				results.add(future);
			} catch (Exception e) {
				System.out.println("Error calculating the prrotein nano-environment for the complex: " + line);
				e.printStackTrace();
			}
		}
		br.close();

		int idx = 0;
		for (Future<ProteinNanoEnvironment> future : results) {
			ProteinNanoEnvironment nanoEnvironment = future.get();
			if(idx++ == 0) {
				String names[] = nanoEnvironment.getDescriptorsNames();
				System.out.print(names[0]);
				for(int i = 1; i < names.length; i++) System.out.print("\t" + names[i]);
			}
			System.out.println();
			Double values[] = nanoEnvironment.getDescriptorsValues();
			System.out.print(Utils.round(values[0],4));
			for(int i = 1; i < values.length; i++) System.out.print("\t" + Utils.round(values[i],4));
		}
	}
}