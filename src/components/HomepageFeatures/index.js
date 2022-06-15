import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';

const FeatureList = [
    {
      title: 'Genomics',
      Svg: require('@site/static/img/undraw_docusaurus_react.svg').default,
      //png:  "https://res.cloudinary.com/dwccfildc/c_limit,f_auto,w_1140/v1637151402/prod/4e093e2e53f12df612635bd5bf54ad58.jpg",
      description: (
        <>

          Learn about genomic technologies and their applications.

        </>
      ),
  },
  {
    title: 'Bioinformatics',
    Svg: require('@site/static/img/undraw_docusaurus_mountain.svg').default,
    description: (
      <>

        A comprehensive and growing set of bioinformatics training resources.

      </>
    ),
  },
  {
    title: 'Statistics',
    Svg: require('@site/static/img/undraw_docusaurus_tree.svg').default,
    description: (
      <>

        Data analysis means statistical analysis. Learn how to model data and interpret it to
        generate robust scientific conclusions.

      </>
    ),
  },
];

function Feature({Svg, title, description}) {
  return (
    <div className={clsx('col col--4')}>
      <div className="text--center">
        <Svg className={styles.featureSvg} role="img" />
      </div>
      <div className="text--center padding-horiz--md">
        <h3>{title}</h3>
        <p>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures() {
  return (
    <section className={styles.features}>
      <div className="container">
        <div className="row">
          {FeatureList.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
