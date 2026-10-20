// Testes isolados de exportação. Não automatizam um navegador nem usam DNA real.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const elements = new Map();
const element = id => {
  if (!elements.has(id)) elements.set(id, {textContent:'',innerHTML:'',hidden:false,value:'',dataset:{},files:[],classList:{toggle(){},add(){},remove(){}},addEventListener(){},click(){},showModal(){},close(){}});
  return elements.get(id);
};
let saved;
const context = vm.createContext({
  document:{getElementById:element,querySelectorAll:()=>[],createElement:()=>({click(){}})},
  window:{addEventListener(){},scrollTo(){}}, location:{hash:''},history:{replaceState(){}},
  URL:{createObjectURL(blob){saved=blob;return 'blob:test';},revokeObjectURL(){}},
  URLSearchParams, Blob:class {constructor(parts,options){this.text=parts.join('');this.type=options.type;}},
  fetch:()=>new Promise(()=>{}),setTimeout:()=>1,clearTimeout(){},setInterval:()=>1,console,
});
vm.runInContext(fs.readFileSync(path.join(__dirname,'../web/app.js'),'utf8'), context);
const synthetic={version:'TEST',summary:{name:'<script>alert(1)</script>',build:'GRCh37',format:'VCF',loci:1,called:1,call_rate:100,sha256:'test',stats:{},warnings:[]},limitations:['Synthetic only'],panel:[{title:'<img src=x onerror=alert(1)>',gene:'TEST',rsid:'rs1',chrom:'1',pos:10,genotype:'A/G',status:'observed',effect_allele:'G',interpretation:'Text & context',evidence:'TEST',source:'javascript:alert(1)',coordinate_source:'https://example.org/'}],clinical:{available:true,metadata:{file_date:'TEST',records:1,source:'https://example.org/',sha256:'test'},total_findings:1,message:'Synthetic',findings:[{gene:'  =1+1',rsid:'rs1',chrom:'1',pos:10,ref:'A',alt:'G',genotype:'A/G',copies:1,classification:'Uncertain_significance',conditions:'<script>bad()</script>',review:'TEST',stars:0,source:'https://example.org/'}]},ancestry:{available:false}};
vm.runInContext(`report=${JSON.stringify(synthetic)};exportHTML();`,context);
assert.ok(saved.text.includes('&lt;script&gt;alert(1)&lt;/script&gt;'));
assert.ok(!saved.text.includes('<script>'));
assert.ok(!saved.text.includes('href="javascript:'));
assert.ok(saved.text.includes("default-src 'none'"));
assert.equal(saved.type,'text/html;charset=utf-8');
vm.runInContext('exportCSV();',context);
assert.ok(saved.text.includes('"\'  =1+1"'));
assert.ok(saved.text.startsWith('\ufeff'));
synthetic.ancestry={available:true,markers:600,chromosomes:22,source:'SYNTHETIC ONLY',method:'Test method',limitation:'No real individual',proportions:[{population:'A',label:'Componente A',percent:60,lower:50,upper:70}],regional_proportions:[{population:'Grupo A',percent:60,lower:50,upper:70}],uncertainty:{replicates:200,blocks:50,method:'Block resampling',limitation:'Conditional on test reference'}};
vm.runInContext(`report=${JSON.stringify(synthetic)};exportHTML();`,context);
assert.ok(saved.text.includes('Faixa de reamostragem'));
assert.ok(saved.text.includes('50,0–70,0%'));
assert.ok(saved.text.includes('200 reamostragens'));
console.log('Exportação HTML: escape, URLs e CSP OK. CSV: proteção de fórmula e UTF-8 OK.');
