package br.embrapa.cnptia.gpbc.plc.descriptors.protein;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.descriptors.contacts.Contact;
import br.embrapa.cnptia.cbi.sdl.descriptors.contacts.ContactsCalculator;
import br.embrapa.cnptia.cbi.sdl.descriptors.prototypes.Descriptors;
import br.embrapa.cnptia.cbi.sdl.utils.Utils;
import br.embrapa.cnptia.gpbc.plc.data.MaxContactsTable;
import br.embrapa.cnptia.gpbc.plc.structure.ProteinLigandComplex;

public class ProteinNanoEnvironment {

	private ProteinLigandComplex proteinLigandComplex;
	private MaxContactsTable maxCon;

	private Descriptors<Contact> contacts;
	private ContactsPerResidue conRes; 

	private List<String> descriptorsNames;
	private List<Double> descriptorsValues;
	private List<IProteinDescriptor> descriptors;

	public ProteinNanoEnvironment(ProteinLigandComplex complex, MaxContactsTable maxCon) throws Exception {
		if(complex == null) 
			throw new NullPointerException("ProteinLigandComplex can not be null!");

		if(maxCon == null) 
			throw new NullPointerException("MaxContactsTable can not be null!");

		this.proteinLigandComplex = complex;
		this.maxCon = maxCon;

		descriptorsNames  = new ArrayList<>();
		descriptorsValues = new ArrayList<>();
		descriptors = new ArrayList<>();

		calculateResidueDescriptors();
	}

	public String[] getDescriptorsNames() {
		String[] array = new String[descriptorsNames.size()];
		array = descriptorsNames.toArray(array);
		return array;
	}

	public Double[] getDescriptorsValues() {
		Double[] array = new Double[descriptorsValues.size()];
		array = descriptorsValues.toArray(array);
		return array;
	}

	private void calculateResidueDescriptors() throws Exception {
		ContactsCalculator conCalc = new ContactsCalculator();
		contacts = conCalc.calculate(proteinLigandComplex.getProtein().getStructure(), null);
		conRes = new ContactsPerResidue(contacts);

		descriptors.add(new DensitySponge(proteinLigandComplex.getProtein()));
		descriptors.add(new EnergyDensity(proteinLigandComplex.getProtein(), conRes));
		descriptors.add(new UnusedContacts(proteinLigandComplex.getProtein(), conRes, maxCon));
		
		Curvature curvature = new Curvature(
				proteinLigandComplex.getProtein(),
				"programs/surfrace/surfrace", 
				"programs/surfrace/radii.txt", 
				1, 1.4,"tmp");
		
		descriptors.add(curvature);
		
		SolventAccessibleSurfaceArea sasa = new SolventAccessibleSurfaceArea(
				proteinLigandComplex.getProtein(), 
				"programs/naccess/naccess", 
				"programs/naccess/vdw.radii",
				1.4, 
				"tmp"
				);
		descriptors.add(sasa);

		for(IProteinDescriptor descriptor : descriptors) 
			descriptorsNames.addAll(Arrays.asList(descriptor.getDescriptorNames()));

	}

	public void calculate(double maxDist) throws Exception {
		Double[] values = new Double[descriptorsNames.size()];
		Arrays.fill(values, 0d);
		descriptorsValues.addAll(Arrays.asList(values));

		for(IResidue residue : proteinLigandComplex.getLigandNeighborResidues()){
			if(!Utils.isAminoAcid(residue.getName())) continue;
			AminoAcid aa = (AminoAcid) residue;
			double dist = proteinLigandComplex.getDistanceFromAAToLigand(aa);
			if(dist <= maxDist) {	
				int idx = 0;
				for(IProteinDescriptor descriptor : descriptors) {
					Double[] aaDescriptorValues = descriptor.getDescriptorValues(aa);
					for(Double aaValue : aaDescriptorValues) {
						Double aux = descriptorsValues.get(idx);
						descriptorsValues.set(idx++, aux + aaValue*Math.exp(-dist));
					}
				}
			}
		}
	}
}
