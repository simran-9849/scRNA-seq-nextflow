// cargo add noodles@0.86.0 --features bam,bgzf,cram,sam
// cargo add noodles-bgzf@0.34.0 --features libdeflate
// cargo add noodles-util@0.56.0 --features alignment

// Original code is from: https://github.com/zaeleus/noodles/discussions/313
use std::{
    fs::File,
    io::{self, BufReader, Read},
    num::NonZero,
    path::Path,
};

use noodles::{
    bam, bgzf, cram,
    sam::{
        self,
        alignment::record::data::field::{Tag, Value},
    },
};

use bstr::{BStr, ByteSlice, BString};
use clap::Parser;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short  = 'i', long = "inputBAM")]
    input_bam: String,

    #[arg(short = 't', long = "threads", default_value_t = 1)]
    threads: usize,
}



fn main() -> io::Result<()> {

    // parse arguments
    let args = Args::parse();

    let worker_count: NonZero<usize> = NonZero::new(args.threads).unwrap();

    let q30 = static_read_as_alignment_record_main(args.input_bam, worker_count)?;

    println!("{q30}");

    Ok(())
}

enum Reader<R> {
    Sam(sam::io::Reader<BufReader<R>>, sam::Record),
    Bam(bam::io::Reader<bgzf::MultithreadedReader<R>>, bam::Record),
    Cram(cram::io::Reader<R>),
}

impl<R> Reader<R>
where
    R: Read + Send + 'static,
{
    fn read_header(&mut self) -> io::Result<sam::Header> {
        match self {
            Reader::Sam(reader, _) => reader.read_header(),
            Reader::Bam(reader, _) => reader.read_header(),
            Reader::Cram(reader) => reader.read_header(),
        }
    }

    fn read_record(&mut self) -> io::Result<Option<&dyn sam::alignment::Record>> {
        match self {
            Reader::Sam(reader, ref mut record) => reader.read_record(record).map(|n| match n {
                0 => None,
                _ => Some(record as &dyn sam::alignment::Record),
            }),
            Reader::Bam(reader, ref mut record) => reader.read_record(record).map(|n| match n {
                0 => None,
                _ => Some(record as &dyn sam::alignment::Record),
            }),
            Reader::Cram(_) => todo!(),
        }
    }
}

impl<R> From<sam::io::Reader<BufReader<R>>> for Reader<R>
where
    R: Read,
{
    fn from(inner: sam::io::Reader<BufReader<R>>) -> Self {
        Self::Sam(inner, sam::Record::default())
    }
}

impl<R> From<bam::io::Reader<bgzf::MultithreadedReader<R>>> for Reader<R> {
    fn from(inner: bam::io::Reader<bgzf::MultithreadedReader<R>>) -> Self {
        Self::Bam(inner, bam::Record::default())
    }
}

impl<R> From<cram::io::Reader<R>> for Reader<R> {
    fn from(inner: cram::io::Reader<R>) -> Self {
        Self::Cram(inner)
    }
}

fn static_read_as_alignment_record_main<P>(src: P, worker_count: NonZero<usize>) -> io::Result<f64>
where
    P: AsRef<Path>,
{
    let mut reader = match src.as_ref().extension().and_then(|ext| ext.to_str()) {
        Some("sam") => File::open(src)
            .map(BufReader::new)
            .map(sam::io::Reader::new)
            .map(Reader::from)?,
        Some("bam") => File::open(src)
            .map(|f| bgzf::MultithreadedReader::with_worker_count(worker_count, f))
            .map(bam::io::Reader::from)
            .map(Reader::from)?,
        _ => unimplemented!(),
    };

    let header = reader.read_header()?;

    let mut q30_vec: Vec<f64> = vec![];

    while let Some(record) = reader.read_record()? {
        if let Ok(Some(q30)) = calculate_record_q30(&header, record) {
            //println!("{}", q30);
            q30_vec.push(q30);
        };

    }
    //println!("q30 length: {}", q30_vec.len());
    if q30_vec.len() > 0 {
        let sum_value = q30_vec.iter().sum::<f64>();
        //println!("sum value: {}", sum_value);
        return Ok(sum_value as f64 / q30_vec.len() as f64)
    } else {
        return Ok(0.0)
    }
}

fn calculate_record_q30(
    _header: &sam::Header,
    record: &dyn sam::alignment::Record,
) -> io::Result<Option<f64>> {

    let flags = record.flags()?;

    // filter records before extracting
    if !flags.is_secondary()
       && !flags.is_supplementary() {

           let linker1 = "GTGA";
           let linker2 = "GACA";

           let mut cr: Option<&BStr> = None;
           let mut ur: Option<&BStr> = None;
           let mut ss: Option<&BStr> = None;
           let mut sq: Option<&BStr> = None;

           for result in record.data().as_ref().iter() {

               match result? {
                   (tag, Value::String(value)) if tag == Tag::new(b'C', b'R') && !value.is_empty() && value.to_string() != "-" => cr = Some(value),
                   (tag, Value::String(value)) if tag == Tag::new(b'U', b'R') && !value.is_empty() && value.to_string() != "-" => ur = Some(value),
                   (tag, Value::String(value)) if tag == Tag::new(b's', b'S') && !value.is_empty() && value.to_string() != "-" => ss = Some(value),
                   (tag, Value::String(value)) if tag == Tag::new(b's', b'Q') && !value.is_empty() && value.to_string() != "-" => sq = Some(value),
                   //(Tag::CELL_BARCODE_ID, Value::String(value)) => {
                   //    cell_barcode_id = Some(value.to_string());
                   //    println!("{:?}", cell_barcode_id);
                   //},
                   _ => {}
               }
               if cr.is_some() && ur.is_some() && ss.is_some() && sq.is_some() {
                   return calc_q30(cr, ur, ss, sq, linker1, linker2)
               }

           };

           // calculation failed
           Ok(None)
       } else {
           // did not pass the filter
           Ok(None)
       }

}


fn calc_q30(
    cr: Option<&BStr>, ur: Option<&BStr>,
    ss: Option<&BStr>, sq: Option<&BStr>,
    linker1: &str,
    linker2: &str
) -> io::Result<Option<f64>> {

        match (cr, ur, ss, sq) {
            (Some(cr), Some(ur), Some(ss), Some(sq)) => {
                // Split cr by underscores
                let parts: Vec<&[u8]> = cr.split_str("_").collect();
                if parts.len() != 3 {
                    return Ok(None);
                }

                let pattern = BString::from([
                    parts[0], linker1.as_bytes(),
                    parts[1], linker2.as_bytes(),
                    parts[2], ur.as_bytes()
                ].concat());

                if let Some(m) = ss.find(&pattern) {
                    let start = m;
                    let end = m + pattern.len();
                    let q30 = &sq[start..end]
                        .as_bytes()
                        .into_iter()
                        .filter(|&x| *x >= 33 + 30 && *x <= 126) // Q >= 30, starts from 33 (illumina 1.8)
                        .count();

                    let percentage = *q30 as f64 / pattern.len() as f64;

                    return Ok(Some(percentage))
                }else{
                    return Ok(None)
                }
            },
            _ => Ok(None)
        }
    }
